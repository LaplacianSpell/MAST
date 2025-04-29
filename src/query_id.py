import os
import pandas as pd
from astroquery.mast import Observations
from astropy.coordinates import SkyCoord
import astropy.units as u
from tqdm import tqdm

# Load the CSV file with object information
# You need to make sure that there is a column called source_id, ra and dec.
script_dir = os.path.dirname(os.path.abspath(__file__))
path_file = os.path.join(script_dir, '../data/metadata/test_query_website.csv')
df = pd.read_csv(path_file)

# Limit to first 3 objects for testing; remove this line for full batch processing
df = df[0:3]

# Store the full JWST observation table for each object
obs_metadata = {}

def query_by_coords(row):
    """
    Query the MAST archive for JWST observations using object coordinates.

    Parameters
    ----------
    row : pandas.Series
        A row from a DataFrame containing at least the following columns:
            - 'ra' : Right Ascension in degrees
            - 'dec' : Declination in degrees
            - 'source_id' : A unique object identifier (used if `criteria_` is 'object')

    Returns
    -------
    pandas.Series
        A Series with the following values:
            - 'Instrument' : The JWST instrument name used in the first matching observation (if any)
            - 'Data_Product_Type' : The type of data product (e.g., 'image', 'spectrum')
            - 'Exp_Time' : Exposure time in seconds
            - 'Status' : 'Public' or 'Private' depending on data availability

        If no JWST observation is found or an error occurs, returns a Series of four None values.
    """
    try:

        ra, dec = row['ra'], row['dec']
        source_ID = row['source_id']
        #create the coordinates using the SkyCoord from astropy library
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')

        #Here we make the query using the coordinates from the previous step
        obs_table = Observations.query_region(coord, radius="0.02 deg")
        #Filter the found observations to select just the JWST observations
        jwst_obs = obs_table[obs_table['obs_collection'] == 'JWST']
        if len(jwst_obs) == 0:
            #If there are not observations, we just fill the columns with None values
            return pd.Series([None]*4)

        # Save the full observation table for later download
        obs_metadata[source_ID] = jwst_obs

        # Extract key metadata fields if available
        # All keywords can be found in: https://mast.stsci.edu/api/v0/_c_a_o_mfields.html
        instrument = jwst_obs['instrument_name'][0] if 'instrument_name' in jwst_obs.colnames else None
        exp_type = jwst_obs['dataproduct_type'][0] if 'dataproduct_type' in jwst_obs.colnames else None
        exp_time = jwst_obs['t_exptime'][0] if 't_exptime' in jwst_obs.colnames else None
        status = 'Public' if jwst_obs['dataRights'][0] == 'PUBLIC' else 'Private'

        #There are 4 new columns for each found object
        return pd.Series([instrument, exp_type, exp_time, status])

    except Exception as e:
        # Return None values if any error occurs during querying
        return pd.Series([None]*4)


def download_jwst_products(obs_table, source_id):
    """
    Download public data products for the first JWST observation in the table.

    Parameters
    ----------o
    obs_table : Table
        An astropy Table returned by Observations.query_*().

    source_id : str
        The source identifier used for naming the download folder.

    Returns
    -------
    None
    """
    try:
        obs_id = obs_table[0]['obsid']
        #Filter for Public observations. You can change the filters in this part.
        if obs_table[0]['dataRights'] != 'PUBLIC':
            print(f"[{source_id}] Observation is private. Skipping download.")
            return

        #This part gets the observation files
        products = Observations.get_product_list(obs_id)
        #if you want ALL data, not just SCIENCE data, chance "SCIENCE" for "ALL"
        manifest = Observations.download_products(products, productType='SCIENCE', mrp_only=False)
        print(f"[{source_id}] Download complete.")
    except Exception as e:
        print(f"[{source_id}] Download failed: {e}")






# Run JWST queries for all sources
tqdm.pandas(desc='Querying JWST...')
df[['Instrument', 'Data_Product_Type', 'Exp_Time', 'Status']] = df.progress_apply(query_by_coords, axis=1)

# Print summary
print('Found Object:', sum(~df['Instrument'].isna()))


## Download products for public observations

#Uncomment the following 4 lines to download the science data for each found sources
#for _, row in df.iterrows():
#    sid = row['source_id']
#    if sid in obs_metadata and row['Status'] == 'Public':
#        download_jwst_products(obs_metadata[sid], sid)
