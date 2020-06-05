# This script will find all the individuals who have a desired pedigree in the .network files from the pedigrees

def searchPedigree(row, outputPath):
    with open("/data100t1/home/grahame/projects/old_data/testALLfiles/EUR_MEGA_ex_Array_BestOfMultipleCalls_Below_Stuttering_A1_updatedrsids_only_noindels_nodups-updated_finalQCfile_cleaned.genome_networks") as pedigreeFile:

        for line in pedigreeFile:

            line = line.split()

            if line[2] == row:
                MyFile = open(
                    outputPath, 'w')

                MyFile.write(line)
                MyFile.write('\n')
                MyFile.close()


def findPedigreeInd(filepath, outputPath):
    with open(filepath) as idFile:

        for row in idFile:
            print(row)
            searchPedigree(row, outputPath)
