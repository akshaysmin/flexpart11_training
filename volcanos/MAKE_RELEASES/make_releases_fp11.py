# Making RELEASES 
# Written by Delia Arnold, adapted by MD Mulder

import os
import numpy as np
import release_parameters as relparams

def makeRELEASES():

    """
    Create RELEASES file
   
    Arguments: parameters: A releases-specific dictionary of parameters   
    """
    if not os.path.isdir(relparams.releaseFileDir ):
        print('makeRELEASES(): directory not found: ' + relparams.releaseFileDir)
        sys.exit()

    latList = [relparams.vlat]
    lonList = [relparams.vlon]
    massList = [1]

    intervals = np.arange(relparams.numRelLevels)
    bottomHeights = relparams.relHeightInt+intervals*relparams.relHeightInt
    topHeights = relparams.relHeightInt+ (intervals+1)*relparams.relHeightInt

    releasesFilename = relparams.releaseFileDir + '/RELEASES'

    with open(releasesFilename, 'w') as R:
        R.write('&RELEASES_CTRL \n')
        R.write(' NSPEC= %3d,\n' % len(relparams.speciesIndices))
        R.write(' SPECNUM_REL= ' +str(relparams.speciesIndices).strip('[]')+',\n')
        R.write ('/\n')
    
        siteCounter = 1 # this is the release counter
        for theLevel, theBottom, theTop in zip(intervals, bottomHeights, topHeights):
            for theSpeciesIdx in relparams.speciesIndices:
                for theLat, theLon, themass in zip(latList, lonList, massList):
                    R.write('&RELEASE\n')
                    R.write(' IDATE1= ' + str(relparams.startTime)[0:8] + ',\n')
                    R.write(' ITIME1= ' + str(relparams.startTime)[8:10] + '0000,\n')
                    R.write(' IDATE2= ' + str(relparams.stopTime)[0:8] + ',\n')
                    R.write(' ITIME2= ' + str(relparams.stopTime)[8:10] + '0000,\n')
                    R.write(' LON1= %9.4f,\n' % theLon)
                    R.write(' LON2= %9.4f,\n' % theLon)
                    R.write(' LAT1= %9.4f,\n' % theLat)
                    R.write(' LAT2= %9.4f,\n' % theLat)
                    R.write(' Z1= %10.3f,\n' % theBottom)
                    R.write(' Z2= %10.3f,\n' % theTop)
                    R.write(' ZKIND= %9d,\n' % relparams.elevCode)
                    string = ' MASS= '
                    for theSpeciesJdx in relparams.speciesIndices:
                        if theSpeciesJdx == theSpeciesIdx:
                            string = string + format(float(themass), '9.4E').replace("+", "")  + ','
                        else:
                            string = string + '0.00,'
                    R.write(string + '\n')
                    R.write(' PARTS= %9d,\n' % relparams.particlesPerRelease)
                    R.write(' COMMENT="ERU_%05d",\n' % siteCounter)
                    R.write('/\n')
                    siteCounter += 1
    
if os.path.isfile(relparams.releaseFileDir+'RELEASES'):
    print('File created: '+relparams.releaseFileDir+'RELEASES')
else:
    print('No file created.')

#----------------------------------------------
makeRELEASES()
