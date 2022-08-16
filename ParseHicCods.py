#!/usr/bin/env python3.6

import os
import re
import argparse
import glob
import gzip

import pandas

from cod.util.log.Logger import Logger
from cod.util.time.Timer import Timer

import numpy

#===============================================================================
DESC_COMMENT = "Script to parse HiC data from Rao et al. 2014 and determine chromatin contacts between co-expressed genes."
SCRIPT_NAME = "ParseHicCods.py"
#===============================================================================

"""
#===============================================================================
@author: Diogo Ribeiro
@date: 03 May 2019
@copyright: Copyright 2019, University of Lausanne
"Script to parse HiC data from Rao et al. 2014 and determine chromatin contacts between co-expressed genes."
#===============================================================================

#===============================================================================
# General plan:
# 1)       
#===============================================================================

#===============================================================================
# Processing notes:
# 1)       
#===============================================================================
"""

class ParseHicCods(object):
    
    ## Set default arguments
    DEFAULT_PARAM_MAPQ_FOLDER = "MAPQG0"
    DEFAULT_PARAM_NORM = ""
    DEFAULT_PARAM_ALREADY_NORM = 0
    DEFAULT_PARAM_BED_FILE = ""
    DEFAULT_PARAM_TOP_INTERACTIONS = 0
    DEFAULT_PARAM_INTERACTION_CUTOFF = 0
    
    # Flag to write normalised matrices to file
    WRITE_NORM_MATRICES = 1
    
    ## Output files 
    # File with normalised HiC contacts   
    OUTPUT_NORM_HIC = "_normalised"   
    OUTPUT_BED_CONTACTS = "_contacts"
        

    def __init__(self, hicFolder, outputFolder, resolution, mapqFolder, \
                     norm, alreadyNorm, bedFile, topInteractions, interactionCutoff, verbosityLevel):

        Timer.get_instance().start_chrono()

        self.hicFolder = hicFolder
        self.outputFolder = outputFolder
        self.resolution = resolution
        self.mapqFolder = mapqFolder
        self.norm = norm
        self.alreadyNorm = alreadyNorm
        self.bedFile = bedFile
        self.topInteractions = topInteractions
        self.interactionCutoff = interactionCutoff

        Logger.get_instance().info(
            "HiC folder: %s\nOutput folder: %s\nResolution: %s\nMAPQ: %s\
            \nNormalisation: %s\nLoad already normalised files: %s\nFolder with CODs: %s\
            \nTop interactions: %s\nInteraction cutoff %s\n" \
              % (self.hicFolder, self.outputFolder, self.resolution, self.mapqFolder,
                 self.norm, self.alreadyNorm, self.bedFile, self.topInteractions, self.interactionCutoff))

        if not os.path.isdir(self.outputFolder) :
            os.mkdir(self.outputFolder)

        Logger.get_instance().set_level( verbosityLevel)


    def read_process_hic_folder(self):
        """
        Read Rao et al. 2014 HiC "*.RAWobserved" for all chromosomes.
        Normalise each result if wanted.
        
        Folder structure:
         - one folder for each chromosome (e.g. chr22)
         - instead that folder there is a folder for each MAPQ
         - inside the MAPQ folder, there is several HiC files (RAWObserved, normalisation files etc)        

        """

        if self.alreadyNorm:
            # Find all already normalised HiC files            
            search = self.hicFolder + "/*.RAWobserved*" + self.norm + "*" 
        else:
            # Find all input HiC files
            search = self.hicFolder + "/chr*/" + self.mapqFolder + "/*.RAWobserved"

        hicFilesToProcess = glob.glob(search)

        Logger.get_instance().info( "read_process_hic_folder: found %s files to process" % (len(hicFilesToProcess)))
        if len( hicFilesToProcess) == 0:
            raise Exception( "read_process_hic_folder: did not find any files to process: %s" % (search))
            
        # Process and normalise matrices
        matrices = {} # key -> chromosome, value -> contact matrix
        for hicFile in sorted(hicFilesToProcess):            
            
            if self.alreadyNorm:                   
                # try to get the chr info from the file name
                fileName = hicFile.split("/")[-1]
                query = re.search( ".*chr([0-9XY]+).*", fileName )
                if query != None:
                    chro = query.group(1)
                    Logger.get_instance().info( "read_process_hic_folder: load chromosome: %s" % (chro))
                    matrices[chro] = self.read_hic_file(hicFile)
                    continue
                else:
                    raise Exception( "read_process_hic_folder: could not automatically retrieve name of chromosome from file: %s" % (hicFile))                    
            
            else:
                # try to get the chr info from the folder path
                chroPath = hicFile.split("/")[-3]
                if "chr" not in chroPath:
                    raise Exception( "read_process_hic_folder: could not automatically retrieve name of chromosome from path: %s" % (hicFile))
                else:
                    chro = chroPath.replace("chr","")
                    Logger.get_instance().info( "read_process_hic_folder: processing chromosome: %s" % (chro))

            if self.norm != ParseHicCods.DEFAULT_PARAM_NORM:
                # Normalise data
                
                # Find normalisation file
                normFile = glob.glob(self.hicFolder + "/chr" + chro + "/" + self.mapqFolder + "/*." + self.norm + "norm" )
                if len(normFile) != 1:
                    raise Exception( "read_process_hic_folder: could not find normalisation file: %s" % (normFile))
                
                # read hic file and normalise
                matrices[chro] = self.read_hic_file(hicFile, self.read_norm_file(normFile[0]))
                
            else:
                # read hic file without normalising
                matrices[chro] = self.read_hic_file(hicFile)

        self.matrices = matrices


    def read_hic_file(self, hicFile, normData = ""):
        """
        Read "*RAWobserved" contact matrix file from Rao et al. 2014.
        File is a sparse matrix, no header, tab-separated.
         
        Example format:
            9400000 9400000 2646.0
            9400000 9450000 71.0
            9450000 9450000 77.0
            9400000 9500000 28.0

        Note: assumes the HiC data is sorted in a way that bin1 coordinate is <= than bin2 coordinate.
        Normally, HiC contacts are only displayed in one sense (bin1 contacts bin2, not bin2 contacts bin1)
        Rao et al. 2014 documentation says the coordinate provided is the left edge of bin.
        """
        
        # Read sparse matrix as dataframe
        if self.alreadyNorm:
            # read already normalised matrix, produced by this script                
            if ".gz" in hicFile: inFile = gzip.open(hicFile, "rt")
            else: inFile = open(hicFile, "r")

            matrix = {} #key -> coord1, value -> dict. key -> coord1, val -> contact            
            for line in inFile:
                line = line.strip()                    
                spl = line.split("\t")
                coord1 = int(spl[0])
                coord2 = int(spl[1])
                if len(spl) > 2:
                    normalisedContact = float(spl[2])
                else:
                    normalisedContact = numpy.nan
                
                if coord1 not in matrix:
                    matrix[coord1] = {}
                matrix[coord1][coord2] = normalisedContact
                
                if self.topInteractions != ParseHicCods.DEFAULT_PARAM_TOP_INTERACTIONS:
                    # top interactions analysis requires storing two-way contacts, i.e. bin1-bin2 and bin2-bin1
                    if coord2 not in matrix:
                        matrix[coord2] = {}
                    matrix[coord2][coord1] = normalisedContact
                
            inFile.close()
                   
        else:
            if ".gz" in hicFile: inFile = gzip.open(hicFile, "rt")
            else: inFile = open(hicFile, "r")

            #read and normalise (if wanted)
            matrix = {} #key -> coord1, value -> dict. key -> coord1, val -> contact            
            
            for line in inFile:
                line = line.strip()                    
                spl = line.split("\t")
                coord1 = int(spl[0])
                coord2 = int(spl[1])
                contact = float(spl[2])

                if len(normData) > 0:            
                    # normalisation according to Rao et al. 2014: 
                    # "To normalize, an entry M_i,j in a *RAWobserved file, divide the entry by the corresponding norm factors for i and j. "
                    # e.g.: 59.0/(1.2988778370674694*1.6080499717941548)       
                    normalisedContact = contact / (normData[coord1]*normData[coord2])
                else:
                    normalisedContact = contact
                
                if coord1 not in matrix:
                    matrix[coord1] = {}
                matrix[coord1][coord2] = normalisedContact

                if self.topInteractions != ParseHicCods.DEFAULT_PARAM_TOP_INTERACTIONS:
                    # top interactions analysis requires storing two-way contacts, i.e. bin1-bin2 and bin2-bin1
                    if coord2 not in matrix:
                        matrix[coord2] = {}
                    matrix[coord2][coord1] = normalisedContact
                
            if ParseHicCods.WRITE_NORM_MATRICES:
                # Write normalised contacts to file
                outFile = open(self.outputFolder + "/" + hicFile.split("/")[-1] + self.norm + ParseHicCods.OUTPUT_NORM_HIC, "w")
                for coord1 in sorted(matrix):
                    for coord2 in sorted(matrix[coord1]):
                        outFile.write("%s\t%s\t%s\n" % ( coord1, coord2, matrix[coord1][coord2]))
                outFile.close()
        
        return matrix
        

    def read_norm_file(self, normFile):
        """
        Read "*norm" normalisation file from Rao et al. 2014.
        File is has a value for each bin, no header, single column.
         
        Example format:
            NaN
            NaN
            0.20942816074662943
            0.03956117423174896
        """

        Logger.get_instance().info( "read_norm_file: normalising %s file" % (normFile))
        
        normValues = open(normFile,"r").readlines()
        normData = {} # key -> bin, value -> normalisation factor            
        bin = 0
        for val in normValues:
            normData[bin] = float(val)
            bin+=self.resolution

        return normData
        

    def read_process_bed_file(self):
        """
        Reads and processes BED file containing pairs of coordinates in which to report HiC contacts.
        Considers only the strand adjusted start position (e.g. TSS of a gene) to determine the relevant HiC bin.
        
        Example format (BED/TSV):
            #chr cisStart cisEnd  cisPheno cisInfo cisStrand corr corrSign rank distance adjustedPval centralPhenotype centralStart centralEnd centralStrand centralInfo
            1       895967  901095  ENSG00000187961.13_4    info   +       0.12268 -       0       762244  0.317   ENSG00000238009.6_6     89295   133723  -       info
            1       134901  139379  ENSG00000237683.5       info       -       0.10008 +       1       5656    0.68    ENSG00000238009.6_6     89295   133723  - info

        The following columns are essential:
            #chr,cisStart,cisEnd,cisStrand,centralStart,centralEnd,centralStrand
            
        Note: chromosome IDs need to match the matrix files (20,22 etc, NOT chr1,chr2)
        Note: assumes that coordinates are given in the "+"/positive strand, i.e. that the start is always equal or lower than the end position.
        """
        
        # Read bed file dataframe
        bedData = pandas.read_csv(self.bedFile, sep="\t")
                
        # Remove certain entries with NA (specific to CODer.py files)
        if "rank" in bedData: 
            if bedData["rank"].dtype != numpy.int64:
                bedData = bedData[bedData["rank"] != "central"]
        bedData = bedData.reset_index(drop=True)

        bins = {} # key -> chromosome, val -> list of bins in the chromosome
        for chro in self.matrices:
            bins[chro] = numpy.array(list(self.matrices[chro].keys()))
        
        missingChro = set() # stores set of missing chromosomes
        contactResults = {} # key -> idx, val -> contacts 

        if self.topInteractions != ParseHicCods.DEFAULT_PARAM_TOP_INTERACTIONS:
            # Results to be added to final BED file
            sameBin = {} # key -> idx, val -> flag whether the pair is in the same bin
            otherBinTop = {} # key -> idx, val -> whether binA is found among top interactors of binB (0, 1 or 2 times/reciprocal)
            thirdBin = {} # key -> idx, val -> count of third bins common to both tops
            thirdBinList = {} # key -> idx, val -> list of the third bins 
            
        for idx,row in bedData.iterrows():
            if idx % 10000 == 0:
                Logger.get_instance().info( "read_process_bed_file: processed %s entries" % (idx))             
                #Timer.get_instance().step( str(idx))
                
            chro,cisStart,cisEnd,cisStrand,centralStart,centralEnd,centralStrand = row[["#chr","cisStart","cisEnd",
                                                                                       "cisStrand","centralStart","centralEnd","centralStrand"]]
            chro = str(chro)
            
            if chro in self.matrices:

                assert cisStart <= cisEnd, "Start position should not be smaller than the end position."
                assert centralStart <= centralEnd,  "Start position should not be smaller than the end position."
                             
                # Get strand-adjusted start position (TSS)
                if cisStrand == "+": coord1 = int(cisStart)
                else: coord1 = int(cisEnd)
                if centralStrand == "+": coord2 = int(centralStart)
                else: coord2 = int(centralEnd)
    
                # sort the coordinates so that we start with the smaller            
                coordFirst, coordSecond = sorted([coord1,coord2])
                            
                coordFirstBin = int(int(coordFirst/self.resolution)*self.resolution)
                coordSecondBin = int(int(coordSecond/self.resolution)*self.resolution)
                
                # Analysis of the top HiC interactions
                if self.topInteractions != ParseHicCods.DEFAULT_PARAM_TOP_INTERACTIONS:
                                        
                    # First check if both genes are in the same bin
                    if coordFirstBin == coordSecondBin:
                        # If both are in the same bin, cannot calculate other metrics
                        sameBin[idx] = 1
                        otherBinTop[idx] = "NA"
                        thirdBin[idx] = "NA"
                        thirdBinList[idx] = "NA"
                    else:
                        sameBin[idx] = 0
                    
                        # check if both bins are present
                        if coordFirstBin not in self.matrices[chro] or coordSecondBin not in self.matrices[chro]:
                            Logger.get_instance().warning( "top interaction analysis: one of the bins not present: %s, %s, %s" % (chro, coordFirstBin, coordSecondBin))     
                            otherBinTop[idx] = "NA"
                            thirdBin[idx] = "NA"
                            thirdBinList[idx] = "NA"
                            contactResults[idx] = "NA"
                            continue
                            
                        # check if there are enough total interactions to analyse in these bins, otherwise do not consider this pair
                        if len(self.matrices[chro][coordFirstBin]) > self.topInteractions and len(self.matrices[chro][coordSecondBin]) > self.topInteractions:

                            if self.topInteractions == -1:
                                # special case, in which we pick all interactions, not the top
                                if self.interactionCutoff == ParseHicCods.DEFAULT_PARAM_INTERACTION_CUTOFF:
                                    coordFirstTop = {key for key, _ in self.matrices[chro][coordFirstBin].items()}
                                    coordSecondTop = {key for key, _ in self.matrices[chro][coordSecondBin].items()}
                                else:
                                    coordFirstTop = {key for key, val in self.matrices[chro][coordFirstBin].items() if val > self.interactionCutoff}
                                    coordSecondTop = {key for key, val in self.matrices[chro][coordSecondBin].items() if val > self.interactionCutoff}                                    
                            else:
                                # Sort HiC interactions for the bin, then store the bin of top X interactions
                                sorted1 = sorted(self.matrices[chro][coordFirstBin].items(), key=lambda kv: kv[1], reverse = True)
                                coordFirstTop = set()
                                for binContact in sorted1:
                                    if len(coordFirstTop) == self.topInteractions: break
                                    if binContact[0] != coordFirstBin: # Exclude interactions to the same bin (self interactions)
                                        coordFirstTop.add(binContact[0])
    
                                sorted2 = sorted(self.matrices[chro][coordSecondBin].items(), key=lambda kv: kv[1], reverse = True)
                                coordSecondTop = set()
                                for binContact in sorted2:
                                    if len(coordSecondTop) == self.topInteractions: break 
                                    if binContact[0] != coordSecondBin: # Exclude interactions to the same bin (self interactions)
                                        coordSecondTop.add(binContact[0])
                                                                    
                                assert len(coordFirstTop) == self.topInteractions
                                assert len(coordSecondTop) == self.topInteractions
                                                            
                            # Analysis 1: check if one of the top interactions is the other bin
                            # count will be 1 if one bin is in the other bin's top, but not vice-versa
                            # count will be 2 if both bins are in each other's top, reciprocally
                            bothInTop = 0
                            for bi in coordFirstTop:
                                if bi == coordSecondBin: bothInTop+=1
                            for bi in coordSecondTop:
                                if bi == coordFirstBin: bothInTop+=1
                            otherBinTop[idx] = bothInTop
    
                            # Analysis 2: check common top interactions to a same third bin
                            binOverlap = coordFirstTop.intersection(coordSecondTop)
                            if len(binOverlap) > 0:
                                thirdBin[idx] = len(binOverlap)
                                thirdBinList[idx] = ",".join([str(i) for i in sorted(binOverlap)])
                            else:
                                thirdBin[idx] = 0
                                thirdBinList[idx] = "NA"
                                                
                        else:
                            Logger.get_instance().warning( "top interaction analysis: not enough HiC \
                            interactions for these bins %s, %s, %s" % (chro, coordFirstBin, coordSecondBin))     
                            otherBinTop[idx] = "NA"
                            thirdBin[idx] = "NA"
                            thirdBinList[idx] = "NA"
                                                    
                ## Results about contacts/interactions                
                if coordFirstBin in self.matrices[chro]:
                    if coordSecondBin in self.matrices[chro][coordFirstBin]:
                        contactResults[idx] = self.matrices[chro][coordFirstBin][coordSecondBin]
                    else:
                        contactResults[idx] = "NA"
                else:
                    contactResults[idx] = "NA"
            else:
                contactResults[idx] = "NA"
                missingChro.add( chro)
                continue

                                
        assert bedData.shape[0] == len(contactResults), "Should have obtain a HiC result for each BED entry"

        Logger.get_instance().info( "read_process_bed_file: processed %s BED entries" % (idx+1))
        Logger.get_instance().info( "read_process_bed_file: chromosomes without HiC data: %s" % (",".join(missingChro)))

        # add contacts to bed dataframe and write to file
        bedData["normalised_contact"] = pandas.Series(contactResults)    
        
        # Write to file
        bedData.to_csv(open(self.outputFolder + "/" + self.bedFile.split("/")[-1] + ParseHicCods.OUTPUT_BED_CONTACTS, "w"), 
                                  sep = "\t", index = False, header = True, float_format='%.3f')

        if self.topInteractions != ParseHicCods.DEFAULT_PARAM_TOP_INTERACTIONS:
            # Write to another file the top interaction results
            bedData["same_bin"] = pandas.Series(sameBin)
            bedData["pair_bin_top"] = pandas.Series(otherBinTop)
            bedData["third_bin"] = pandas.Series(thirdBin)
            bedData["third_bin_list"] = pandas.Series(thirdBinList)
            if self.interactionCutoff == ParseHicCods.DEFAULT_PARAM_INTERACTION_CUTOFF:
                outputFileName = self.outputFolder + "/" + self.bedFile.split("/")[-1] + ParseHicCods.OUTPUT_BED_CONTACTS + "_top" + str(self.topInteractions)
            else:
                outputFileName = self.outputFolder + "/" + self.bedFile.split("/")[-1] + ParseHicCods.OUTPUT_BED_CONTACTS + "_top" + str(self.topInteractions) + "_cutoff" + str(self.interactionCutoff)
                
            bedData.to_csv(open(outputFileName, "w"), sep = "\t", index = False, header = True, float_format='%.3f')
            
        return bedData


    def run(self):
        """
        Run functions in order
        """

        Timer.get_instance().step( "Reading matrix file.." )        
        self.read_process_hic_folder()

        if self.bedFile != ParseHicCods.DEFAULT_PARAM_BED_FILE:
            Timer.get_instance().step( "Reading bed file.." )        
            self.read_process_bed_file()


if __name__ == "__main__":

    try:
    
        # Start chrono
        print ("STARTING " + SCRIPT_NAME)

        #===============================================================================
        # Get input arguments
        #===============================================================================
        parser = argparse.ArgumentParser(description= DESC_COMMENT) 
    
        # positional args
        parser.add_argument('hicFolder', metavar='hicFolder', type=str,
                             help='Folder with HiC data from Rao et al. 2014. Choose a resolution and provide the folder containing all chromosome subfolders\
                              (e.g. "GM12878_combined/5kb_resolution_intrachromosomal".')
        parser.add_argument('outputFolder', metavar='outputFolder', type=str,
                             help='Folder where all output files will be written.')
        parser.add_argument('resolution', metavar='resolution', type=int,
                             help='Resolution of HiC data in base-pairs (bp). This has to match the provided folder. E.g. 5000.')
        parser.add_argument('--mapqFolder', metavar='mapqFolder', type=str, default = ParseHicCods.DEFAULT_PARAM_MAPQ_FOLDER,
                             help='Choice of available HiC data MAPQ, based on the subfolder name (e.g. "MAPQG0" or "MAPQGE30") (default: MAPQG0).')        
        parser.add_argument('--norm', metavar='norm', type=str, default = ParseHicCods.DEFAULT_PARAM_NORM, 
                             help='Choice of the normalisation procedure, based on files present. Choices: "KR", "VC", "SQRTVC" (default: OFF).')
        parser.add_argument('--alreadyNorm', metavar='alreadyNorm', type=str, default = ParseHicCods.DEFAULT_PARAM_ALREADY_NORM, 
                             help='Whether to load already normalised HiC data. Use together with --norm to pick the right files. (default: OFF).')
        parser.add_argument('--bedFile', metavar='bedFile', type=str, default = ParseHicCods.DEFAULT_PARAM_BED_FILE, 
                             help='Optional BED file containing coordinates in which contacts (bin) will be reported. \
                             File should contain the following columns (in header): #chr,cisStart,cisEnd,cisStrand,centralStart,centralEnd,centralStrand.\
                             Assumes both pairs of coordinates are in the same chromosome.')
        parser.add_argument('--topInteractions', metavar='topInteractions', type=int, default = ParseHicCods.DEFAULT_PARAM_TOP_INTERACTIONS, 
                             help='Output results of picking X top HiC contacts. If -1, pick all HiC contacts. (default: 0/OFF).')
        parser.add_argument('--interactionCutoff', metavar='interactionCutoff', type=float, default = ParseHicCods.DEFAULT_PARAM_INTERACTION_CUTOFF, 
                             help='HiC interaction cutoff. Works together with --topInteractions == -1, picking all interactions above given value. (default: 0/OFF).')
        parser.add_argument('--verbosityLevel', metavar='verbosityLevel', type=str, default = "debug", 
                             choices = ["debug", "info", "warning", "error", "critical", "fatal"],
                             help='Level of verbosity. Choices: "debug", "info", "warning", "error", "critical", "fatal"')
                
        #gets the arguments
        args = parser.parse_args( ) 
    
        #===============================================================================
        # Run analysis / processing
        #===============================================================================

        # Initialise class    
        run = ParseHicCods( args.hicFolder, args.outputFolder, args.resolution, args.mapqFolder, \
                     args.norm, args.alreadyNorm, args.bedFile, args.topInteractions, args.interactionCutoff, args.verbosityLevel)

        run.run()

        # Stop the chrono      
        Timer.get_instance().stop_chrono( "FINISHED " + SCRIPT_NAME )

    # Use RainetException to catch errors
    except Exception as e:
        Logger.get_instance().error( "Error during execution of %s. Aborting :\n" % SCRIPT_NAME + str(e))

