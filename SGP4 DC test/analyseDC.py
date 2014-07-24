# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 15:03:09 2014

Parse output files from testdc.cpp and create plots that enable accuracy and 
efficiency of the SGP4 Differential Correction (Orbit Determination) to be analysed.
@author: al11g09@soton.ac.uk
@version: 1.0.0
@since 23-07/2013
"""
import numpy, matplotlib.pyplot, matplotlib.pylab, matplotlib

""" Parse the data from the second output file that contains comparison between the initial and final TLEs' state vectors. """
with open("testdc2.out", "r") as InFile2:
    inFile2Lines = InFile2.readlines()
    
inFile2Lines = inFile2Lines[8:] # Get rid of the header.

initialDifferences = [] # Differences in position at the beginning of iterations.
finalDifferences = [] # Differences in position at the end of iterations.
iterations = [] # Number of iterations.
for line in inFile2Lines: # Extract the numbers from every line.
    try:
        temp = [float(s) for s in line.split()] # Lis tof numbers: iter   norad      inti diff (m)  final diff    a           e          i     apAlt   prAlt magnitudeOfCovarianceMatrix
        initialDifferences.append( temp[2] )
        finalDifferences.append( temp[3] )
        iterations.append( temp[0] )
    except ValueError:
        print "Value error for: "+line
    
initialDifferences = numpy.array( initialDifferences ) # Change type to numpy.Array for ease of handling.
finalDifferences = numpy.array( finalDifferences )
initialDifferences.sort()
finalDifferences.sort()
IterationsCapReached = iterations.count( max(iterations) ) # Number of times the imposed limit of iterations has been reached.
iterations = numpy.array( iterations )

" Reference data from D. Vallado's paper. "
RefCounts = numpy.array([16198, 16198, 4, 2, 6, 20, 5], dtype=numpy.float64 )
RefCounts = RefCounts/(RefCounts.sum()-RefCounts[0]) # Normalise.
RefBins = numpy.array([1e-6, 1, 10, 100, 1000, 10000, 1e5], dtype=numpy.float64)

" Find the cases where the iterations diverged, i.e. when the final difference was larger than the original. "
findArray = initialDifferences - finalDifferences # Entries of findArray will be -ve if the iterations diverged.
itemsIndices = numpy.where( findArray<0.0 ) # See whee the iterations diverged.

""" Plot the distribution of TLE reconstruction accuracies. """
ticksFontSize = 12
labelsFontSize = 16
titleFontSize = 20

matplotlib.rc('xtick', labelsize=ticksFontSize) 
matplotlib.rc('ytick', labelsize=ticksFontSize)

fig = matplotlib.pyplot.figure()
fig.suptitle(r'$Precision\ of\ TLE\ reproduction\ algorithm\ by\ David\ Vallado$', fontsize=titleFontSize)         
ax = fig.gca()
# Add extra explanatory text.
matplotlib.pyplot.text(0.5, 1.04, r'${:.4f}\ average\ iterations,\ number\ of\ times\ iterations\ limit\ was\ reached\ was\ {}$'.format(iterations.mean(),IterationsCapReached),
         horizontalalignment='center',
         fontsize=labelsFontSize,
         transform = ax.transAxes)
matplotlib.pyplot.text(0.5, 1.005, r'$Number\ of\ times\ iterations\ diverged\ was\ {}$'.format(len(itemsIndices[0])),
         horizontalalignment='center',
         fontsize=labelsFontSize,
         transform = ax.transAxes)
         
ax.plot(RefBins, RefCounts, 'r-', label=r'$Reference$')   

ax.yaxis.set_ticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1]) # This will be a normalised histogram so can set limits like this - data will bever exceed 1.0.

ax.set_ylim([0, 1.1])
ax.set_xscale('log')
ax.grid(b=True, which='major', color='k', linestyle='-')
ax.grid(b=True, which='minor', color='gray', linestyle='--')

n, bins, patches = ax.hist(finalDifferences, 1000, normed=1, histtype='step', cumulative=True, label=r'$Implementation$')

ax.set_xlabel(r'$Difference\ between\ original\ and\ reproduced\ TLE\ (m)$',fontsize=labelsFontSize)
ax.set_ylabel(r'$Fraction\ of\ the\ total\ sample$',fontsize=labelsFontSize)
matplotlib.pyplot.legend()

""" Plot the number of iterations executed to reproduce the TLEs. """
fig2 = matplotlib.pyplot.figure()
fig2.suptitle(r'$Histogram\ of\ the\ number\ of\ iterations\ done\ do\ reproduce\ the\ TLEs$', fontsize=titleFontSize)         
ax2 = fig2.gca()

ax2.grid(b=True, which='major', color='k', linestyle='-')
ax2.grid(b=True, which='minor', color='gray', linestyle='--')

n, bins, patches = ax2.hist(iterations, int(iterations.max()), normed=1, histtype='bar', cumulative=False, label=r'$Number\ of\ iterations$')

ax2.set_xlabel(r'$Number\ of\ iterations$',fontsize=labelsFontSize)
ax2.set_ylabel(r'$Fraction\ of\ a\ sample\ of\ {}\ objects$'.format( iterations.size ),fontsize=labelsFontSize)
matplotlib.pyplot.legend()

fig2.show()
