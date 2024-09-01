import time
start_time = time.time()
import argparse
import io
import logging
import re
import warnings
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
from math import log10

import pandas as pd
import requests
from pyasassn.client import SkyPatrolClient

# Mute UserWarnings
warnings.filterwarnings("ignore", category=UserWarning)
logging.basicConfig(level=logging.INFO)


class BaseDataLoader(ABC):
    """
    Abstract base class for data loaders.

    Attributes
    ----------
    ra : float
        Right ascension of the source.
    dec : float
        Declination of the source.
    c : float
        Speed of light in cm/s.
    Rv : float
        Reddening parameter.
    lambda_and_const_values : dict
        Dictionary of filter names and their corresponding wavelength and constant values.
    filter_names : dict
        Dictionary mapping filter IDs to filter names.
    """

    # Constants
    c = 2.9979e10  # Speed of light in cm/s
    Rv = 3.1
    lambda_and_const_values = {
        "psy": (9620, log10(3631) - 23),
        "psz": (8665, log10(3631) - 23),
        "psi": (7522, log10(3631) - 23),
        "psr": (6176, log10(3631) - 23),
        "psg": (4810, log10(3631) - 23),
        'psy_ztf': (9620, log10(3631) - 23),
        'psz_ztf': (8660, log10(3631) - 23),
        'psi_ztf': (7520, log10(3631) - 23),
        'psr_ztf': (6170, log10(3631) - 23),
        'psg_ztf': (4810, log10(3631) - 23),
        "V": (5500, log10(3640) - 23),
        "g": (4810, log10(3631) - 23),
    }
    filter_names = {1: 'psg', 2: 'psr', 3: 'psi', 4: 'psz', 5: 'psy'}

    def __init__(self, ra, dec):
        """
        Initialize the data loader.

        Parameters
        ----------
        ra : float
            Right ascension of the source.
        dec : float
            Declination of the source.
        """
        self.ra = ra
        self.dec = dec

    @abstractmethod
    def query_data(self):
        """Query the data source."""
        pass

    @abstractmethod
    def process_data(self, data):
        """Process the data."""
        pass

    @staticmethod
    def format_to_3e(x):
        """
        Format a number to 3 decimal scientific notation.

        Parameters
        ----------
        x : float
            The number to format.

        Returns
        -------
        str
            The number formatted in scientific notation with 3 decimal places.
        """
        return f"{x:.3e}"

    @staticmethod
    def calc_nh_value_requests(ra, dec):
        """
        Calculate NH value using requests to an external service.

        Parameters
        ----------
        ra : float
            Right ascension of the source.
        dec : float
            Declination of the source.

        Returns
        -------
        float
            The NH value retrieved from the service.
        """
        entry = f"{ra},{dec}"
        url = "https://heasarc.gsfc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl"

        params = {
            'Entry': entry,
            'NR': 'GRB/SIMBAD+Sesame/NED',
            'CoordSys': 'Equatorial',
            'equinox': '2000',
            'radius': '0.1',
            'usemap': '0'
        }

        # Make a GET request with the specified parameters
        response = requests.get(url, params=params)
        # Check the response status code
        if response.status_code == 200:
            pattern = r"Weighted average nH \(cm\*\*-2\)\s+(\d+\.\d+E[+-]\d+)"
            match = re.search(pattern, response.text)

            if match:
                nh_value = match.group(1)
                return float(nh_value)
            else:
                print("nH value not found in the response.")
                return None
        else:
            print(f"Request failed with status code: {response.status_code}")
            return None

    def find_lambda_and_const(self, filter_name):
        """
        Find lambda and constant values for a given filter.

        Parameters
        ----------
        filter_name : str
            The name of the filter for which lambda and constant values are to be found.

        Returns
        -------
        lambda_ : float
            The lambda value for the given filter.
        const_ : float
            The constant value for the given filter.

        Raises
        ------
        ValueError
            If the filter is not supported, a ValueError is raised with the message "mag2flux: Filter not supported".
        """
        if filter_name not in self.lambda_and_const_values.keys():
            print("mag2flux: Filter not supported")
            raise ValueError("mag2flux: Filter not supported")
        else:
            lambda_, const_ = self.lambda_and_const_values.get(filter_name, (None, None))
        if not lambda_ and not const_:
            print("mag2flux: Filter not supported")
            return
        else:
            return lambda_, const_

    def calculate_extinction(self, lambda_, nh_value):
        """
        Calculates the extinction value based on the input wavelength and NH value.

        Parameters
        ----------
        lambda_ : float
            Wavelength in Angstroms.
        nh_value : float
            NH value in cm^-2.

        Returns
        -------
        a_band: float
            The calculated extinction value.
        """
        x = (10000 / lambda_)
        av = max(self.Rv * (-0.055 + nh_value * 1.987e-22), 0)
        ebv = av / self.Rv

        if 1.1 > x >= 0.3:
            a_band = (0.574 * (x ** 1.61) - 0.527 * (x ** 1.61) / self.Rv) * av
        elif 3.3 > x >= 1.1:
            y = x - 1.82
            aa = 1 + (0.17699 * y) - (0.50447 * y ** 2) - (0.02427 * y ** 3) + (0.72085 * y ** 4) + (
                    0.01979 * y ** 5) - (
                         0.77530 * y ** 6) + (0.32999 * y ** 7)
            bb = 1.41338 * y + (2.28305 * y ** 2) + (1.07233 * y ** 3) - (5.38434 * y ** 4) - (0.62251 * y ** 5) + (
                    5.30260 * y ** 6) - (2.09002 * y ** 7)
            a_band = (aa + (bb / self.Rv)) * av
        elif 10.0 >= x >= 3.3:
            c2 = -0.824 + (4.717 / self.Rv)
            c1 = 2.03 - (3.007 * c2)
            dx = (x * x) / (((x ** 2) - (4.596 ** 2)) ** 2 + ((x * 0.99) ** 2))
            px = (0.5392 * ((x - 5.9) ** 2)) + (0.05644 * ((x - 5.9) ** 3))
            if x <= 5.9:
                px = 0
            a_band = (c1 + c2 * x + 3.23 * dx + 0.41 * px) * ebv + av
        else:
            a_band = 0

        return a_band

    def calculate_flux_and_frequency(self, nh, filter_name, magnitude=0, flux=None, flux_err=None):
        lambda_, const = self.find_lambda_and_const(filter_name)
        frequency = self.c / (lambda_ * 1.e-8)
        a_band = self.calculate_extinction(lambda_, nh)
        extinction_value = 10. ** const / 10. ** (0.4 * a_band + const)

        if flux is not None and flux_err is not None:
            flux = flux * 1.e-23 * frequency / extinction_value
            flux_err = flux_err * 1.e-23 * frequency
            return pd.Series([flux, flux_err, frequency])
        else:
            if magnitude == 0:
                flux = 0
            else:
                flux = 10. ** (-0.4 * (magnitude - a_band) + const) * frequency
            return pd.Series([flux, frequency])

    def add_flux_and_frequency(self, data, nh_value, filter_column='phot_filter', mag_column='mag',
                               mag_err_column='mag_err', psf_flux_column='psfFlux', psf_flux_err_column='psfFluxErr'):
        if psf_flux_column in data.columns and psf_flux_err_column in data.columns:
            data = pd.concat([data, data.apply(
                lambda x: self.calculate_flux_and_frequency(nh_value, x[filter_column], flux=x[psf_flux_column],
                                                            flux_err=x[psf_flux_err_column]), axis=1).rename(
                columns={0: "Flux", 1: "Flux_err", 2: "Frequency"})], axis=1)
        else:
            data = pd.concat([data, data.apply(
                lambda x: self.calculate_flux_and_frequency(nh_value, x[filter_column], magnitude=x[mag_column]),
                axis=1).rename(
                columns={0: "Flux", 1: "Frequency"})], axis=1)
            data = pd.concat([data, data.apply(
                lambda x: self.calculate_flux_and_frequency(nh_value, x[filter_column],
                                                            magnitude=x[mag_column] + x[mag_err_column]),
                axis=1).rename(columns={0: "Flux_err_minus", 1: "drop_column"})], axis=1)
            data["Flux_err"] = data.Flux - data.Flux_err_minus
        return data

    def format_data_columns(self, data):
        data[['Frequency', 'Flux', 'Flux_err']] = data[['Frequency', 'Flux', 'Flux_err']].apply(
            lambda col: col.map(self.format_to_3e))
        return data

    @staticmethod
    def finalize_dataframe(data, catalog_name, reference):
        data['Flag'] = ''
        data['Catalog'] = catalog_name
        data['Reference'] = reference
        data = data.loc[:, ["Frequency", "Flux", "Flux_err", "MJD_start", "MJD_end", "Flag", "Catalog", "Reference"]]
        data.columns = ["freq. ", "flux ", "err_flux ", "MJD_start ", "MJD_end ", "flag", "catalog", "reference"]
        return data

    def query_process_save(self):
        data = self.query_data()
        if not data.empty:
            processed_data = self.process_data(data)
        else:
            processed_data = pd.DataFrame(
                columns=["freq.", "flux", "err_flux", "MJD_start", "MJD_end", "flag", "catalog", "reference"])
        return processed_data


class ASASSNDataLoader(BaseDataLoader):
    """
    Data loader for ASAS-SN data.
    """

    def __init__(self, ra, dec, radius):
        """
        Initialize the ASASSNDataLoader with the given coordinates and search radius.

        Parameters
        ----------
        ra : float
            Right ascension of the source.
        dec : float
            Declination of the source.
        radius : float
            Search radius for the query.
        """
        super().__init__(ra, dec)
        self.radius = radius
        self.catalog_name = 'ASAS-SN-LC'
        self.reference = 'Hart, K., Shappee, B. J., et al. 2023, arXiv:2304.03791'

    def query_data(self):
        """
        Query ASAS-SN data for the given coordinates and radius.

        Returns
        -------
        pd.DataFrame
            Data containing the queried data.
        """
        client = SkyPatrolClient(verbose=False)
        try:
            # Perform a cone search using the SkyPatrolClient
            data = client.cone_search(self.ra, self.dec, radius=self.radius, download=True, threads=10,
                                      units="arcsec").data
        except ValueError:
            logging.warning(
                f'There are no available ASAS-SN observations with Ra={self.ra} Dec={self.dec} coordinates.')
            # Create an empty DataFrame if no data is available
            data = self.create_empty_dataframe()
        return data

    @staticmethod
    def create_empty_dataframe():
        """
        Create an empty DataFrame with the required columns.

        Returns
        -------
        pd.DataFrame
            An empty DataFrame with predefined columns.
        """
        columns = ['asas_sn_id', 'jd', 'flux', 'flux_err', 'mag', 'mag_err', 'limit', 'fwhm', 'image_id', 'camera',
                   'quality', 'phot_filter']
        return pd.DataFrame(columns=columns)

    def process_data(self, data):
        """
        Process the ASAS-SN data.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Processed DataFrame.
        """
        # Filter data based on quality and magnitude error
        data = self.filter_data(data)
        # Add MJD columns
        data = self.add_mjd_columns(data)
        # Calculate NH value
        nh_value = self.calc_nh_value_requests(self.ra, self.dec)
        # Add flux and frequency columns
        data = self.add_flux_and_frequency(data, nh_value)
        # Format data columns
        data = self.format_data_columns(data)
        # Finalize the DataFrame
        data = self.finalize_dataframe(data, self.catalog_name, self.reference)
        logging.info("ASAS-SN data formatting completed.")
        return data

    @staticmethod
    def filter_data(data):
        """
        Filter the ASAS-SN data based on quality and magnitude error.

        Parameters
        ----------
        data : pd.DataFrame
            DataFrame containing the queried data.

        Returns
        -------
        pd.DataFrame
            Filtered DataFrame.
        """
        return data[(data["quality"] == "G") & (data["mag_err"] < 20)]

    @staticmethod
    def add_mjd_columns(data):
        """
        Add MJD start and end columns to the data.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Data with MJD start and end columns.
        """
        data = data.copy()  # Make a copy to avoid SettingWithCopyWarning
        data.loc[:, "MJD_start"] = data["jd"] - 2400000.5
        data.loc[:, "MJD_start"] = data["MJD_start"].round(6)
        data.loc[:, "MJD_end"] = data["MJD_start"]
        return data


class PanSTARRDataLoader(BaseDataLoader):
    """
    Data loader for Pan-STARRS data.
    """

    def __init__(self, ra, dec, radius):
        """
        Initialize the PanSTARRDataLoader with the given coordinates and search radius.

        Parameters
        ----------
        ra : float
            Right ascension of the source.
        dec : float
            Declination of the source.
        radius : float
            Search radius for the query.
        """
        super().__init__(ra, dec)
        self.radius = radius  # Default radius, can be adjusted as needed
        self.catalog_name = "Pan-STARRS-LC"
        self.reference = "Chambers, K. C., Magnier, E. A., et al. 2016, arXiv:1612.05560"

    def query_data(self):
        """
        Query Pan-STARRS data for the given coordinates and radius.

        Returns
        -------
        pd.DataFrame
            Data containing the queried data.
        """
        base_path = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/detection.csv"
        columns_to_get = "columns=[ObjName,obsTime,filterID,psfFlux,psfFluxErr]"
        url = f"{base_path}?ra={self.ra}&dec={self.dec}&radius={self.radius / 3600}&pagesize=500001&{columns_to_get}"
        response = requests.get(url)
        if response.status_code == 200:
            # Read the response content into a DataFrame
            data = pd.read_csv(io.StringIO(response.content.decode('utf-8')))
            if data.empty:
                logging.warning(
                    f'There are no available Pan-STARR observations with Ra={self.ra} Dec={self.dec} coordinates.')
                # Create an empty DataFrame if no data is available
                data = self.create_empty_dataframe()
        else:
            logging.warning(f'Failed to fetch Pan-STARR data for RA: {self.ra}, DEC: {self.dec}')
            # Create an empty DataFrame if the request fails
            data = self.create_empty_dataframe()
        return data

    @staticmethod
    def create_empty_dataframe():
        """
        Create an empty DataFrame with the required columns.

        Returns
        -------
        pd.DataFrame
            An empty DataFrame with predefined columns.
        """
        columns = ['ObjName', 'obsTime', 'filterID', 'psfFlux', 'psfFluxErr']
        return pd.DataFrame(columns=columns)

    def process_data(self, data):
        """
        Process the Pan-STARRS data.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Processed DataFrame.
        """
        # Filter data based on PSF flux
        data = self.filter_data(data)
        # Add MJD columns
        data = self.add_mjd_columns(data)
        # Calculate NH value
        nh_value = self.calc_nh_value_requests(self.ra, self.dec)
        # Add filter name based on filter ID
        data['filter_name'] = data['filterID'].apply(lambda x: self.filter_names[x])
        # Add flux and frequency columns
        data = self.add_flux_and_frequency(data, nh_value, filter_column='filter_name', mag_column='psfFlux',
                                           mag_err_column='psfFluxErr', psf_flux_column='psfFlux',
                                           psf_flux_err_column='psfFluxErr')
        # Format data columns
        data = self.format_data_columns(data)
        # Finalize the DataFrame
        data = self.finalize_dataframe(data, self.catalog_name, self.reference)
        logging.info("Pan-STARRS data formatting completed.")
        return data

    @staticmethod
    def filter_data(data):
        """
        Filter the Pan-STARRS data based on PSF flux.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Filtered DataFrame.
        """
        return data[data['psfFlux'] > 0]

    @staticmethod
    def add_mjd_columns(data):
        """
        Add MJD start and end columns to the data.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Data with MJD start and end columns.
        """
        data = data.copy()  # Make a copy to avoid SettingWithCopyWarning
        data.loc[:, "MJD_start"] = data["obsTime"]
        data.loc[:, "MJD_start"] = data["MJD_start"].round(6)
        data.loc[:, "MJD_end"] = data["MJD_start"]
        return data


class ZTFDataLoader(BaseDataLoader):
    """
    Data loader for ZTF data.
    """

    def __init__(self, ra, dec, radius):
        """
        Initialize the ASASSNDataLoader with the given coordinates and search radius.

        Parameters
        ----------
        ra : float
            Right ascension of the source.
        dec : float
            Declination of the source.
        radius : float
            Search radius for the query.
        """
        super().__init__(ra, dec)
        self.radius = radius
        self.catalog_name = "ZTF-LC"
        self.reference = "Bellm et al. 2019 PASP 131 018002"

    def query_data(self):
        """
        Query ZTF data for the given coordinates and radius.

        Returns
        -------
        pd.DataFrame
            Data containing the queried data.
        """
        prefix_of_url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"
        #0.0014 -> 5 arcsecs radius
        url = f"{prefix_of_url}POS=CIRCLE%20{self.ra}%20{self.dec}%200.0014&BANDNAME=g,r,i&FORMAT=csv"
        response = requests.get(url)
        if response.status_code == 200:
            data = pd.read_csv(io.StringIO(response.content.decode('utf-8')))
            data = data.loc[:, ['mjd', 'filtercode', 'mag', 'magerr']]
            if data.empty:
                logging.warning(
                    f'There are no available ZTF observations with Ra={self.ra} Dec={self.dec} coordinates.')
                data = self.create_empty_dataframe()
        else:
            logging.warning(f'Failed to fetch ZTF data for RA: {self.ra}, DEC: {self.dec}')
            data = self.create_empty_dataframe()
        return data

    @staticmethod
    def create_empty_dataframe():
        """
        Create an empty DataFrame with the required columns.

        Returns
        -------
        pd.DataFrame
            An empty DataFrame with predefined columns.
        """
        columns = ['mjd', 'filtercode', 'mag', 'magerr']
        return pd.DataFrame(columns=columns)

    def process_data(self, data):
        """
        Process the ZTF data.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Processed DataFrame.
        """
        # Correct filter names
        data = self.correct_filter_name(data)
        # Add MJD columns
        data = self.add_mjd_columns(data)
        # Calculate NH value
        nh_value = self.calc_nh_value_requests(self.ra, self.dec)
        # Add flux and frequency columns
        data = self.add_flux_and_frequency(data, nh_value, filter_column='filtercode', mag_column='mag',
                                           mag_err_column='magerr')
        # Format data columns
        data = self.format_data_columns(data)
        # Finalize the DataFrame
        data = self.finalize_dataframe(data, self.catalog_name, self.reference)
        logging.info("ZTF data formatting completed.")
        return data

    @staticmethod
    def correct_filter_name(data):
        """
        Correct the filter names in the data.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Data with corrected filter names.
        """
        data['filtercode'] = "ps" + data['filtercode'].str[1] + "_ztf"
        return data

    @staticmethod
    def add_mjd_columns(data):
        """
        Add MJD start and end columns to the data.

        Parameters
        ----------
        data : pd.DataFrame
            Data containing the queried data.

        Returns
        -------
        pd.DataFrame
            Data with MJD start and end columns.
        """
        data = data.copy()  # Make a copy to avoid SettingWithCopyWarning
        data.loc[:, "MJD_start"] = data["mjd"].round(6)
        data.loc[:, "MJD_end"] = data["MJD_start"]
        return data


def run_all_data_loaders(ra, dec):
    """
    Run all data loaders concurrently.

    Parameters
    ----------
    ra : float
        Right ascension of the source.
    dec : float
        Declination of the source.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame containing results from all data loaders.
    """
    loaders = [ASASSNDataLoader(ra, dec, 5), PanSTARRDataLoader(ra, dec, 5), ZTFDataLoader(ra, dec, 5)]
    results = []

    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(loader.query_process_save) for loader in loaders]
        for future in as_completed(futures):
            try:
                result = future.result()
                if not result.empty:
                    results.append(result)
            except Exception as e:
                print(f"Error processing source: {e}")

    return pd.concat(results, ignore_index=True)


def main():
    """
    Main function to parse arguments and run the appropriate data loader.

    Returns
    -------
    pd.DataFrame
       The resulting data from the selected data loader.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--ra', type=float, help='Right ascension of the source (in degrees).', required=True)
    parser.add_argument('--dec', type=float, help='Declination of the source (in degrees).', required=True)
    parser.add_argument('--source', type=str, help='Source of data (asas_sn, pan_starr, ztf, all).', required=True)

    args = parser.parse_args()
    ra = args.ra
    dec = args.dec
    source = args.source

    if source == 'all':
        data = run_all_data_loaders(ra, dec)
        data.to_csv('optical_data.csv')
    else:
        if source.lower() == 'asas_sn':
            data_loader = ASASSNDataLoader(ra, dec, radius=3)
        elif source.lower() == 'pan_starr':
            data_loader = PanSTARRDataLoader(ra, dec, radius=3)
        elif source.lower() == 'ztf':
            data_loader = ZTFDataLoader(ra, dec, radius=3)
        else:
            print("Invalid source specified.")
            return
        data = data_loader.query_process_save()
        data.to_csv('optical_data.csv')
        #data.to_csv(f'{source}.csv', index=False)
    return data


if __name__ == '__main__':
    main()
    print(f'It took: {time.time()-start_time} seconds ')
