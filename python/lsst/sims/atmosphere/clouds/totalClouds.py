import numpy

class TotalCloudExtinction():
    '''This class is used to generate total cloud extinction for a set of pointings, using appropriate
    long-term weather trends to generate this extinction.
    Based on data from XXXXX '''

    def __init__(self):
        '''Instantiate the TotalCloudExtinction object.'''
        self.total_extinction = None
        return

    def generateTotalClouds(self, opsim_visits):
        '''Generate the atmospheric parameters over time.
        
        The input parameter, 'opsim_visits', is a dictionary of numpy arrays containing opsim visit information
        as opsim_visits['expmjd'] (a numpy array of expmjd's), 'ra', 'dec', and 'obshistid' (at a minimum).
        These arrays are already sorted order of increasing expmjd. 
        Returns a list/numpy array with the total cloud extinction value for each opsim visit (in the same order).'''  
        # Take the data from Pachon or wherever else we can get a cloud extinction curve, fit it, and then
        # use that (with some randomness) to generate the cloud extinction for these opsim visits.
        # Must take into account the fact that opsim is already including some cloudiness which results in missed
        # visits. How to sync up opsim weather and this cloud extinction pattern??

        return self.total_extinction
    
