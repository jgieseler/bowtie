#!/usr/bin/env python3
"""
This file contains the Spectra class that is used to run bowtie analysis.
"""
__author__ = "Christian Palmroos"
__credits__ = ["Christian Palmroos", "Philipp Oleynik"]

import numpy as np

from . import bowtie_util
from . import bowtie_calc 
from . import validations as validate

class Spectra:
    """
    Contains the information about what kind of spectra are considered for 
    the bowtie analysis.
    """

    def __init__(self, gamma_min:float, gamma_max:float, gamma_steps:int=100,
                 cutoff_energy:float=0.0) -> None:

        self.gamma_min: float = gamma_min
        self.gamma_max: float = gamma_max
        self.gamma_steps: int = gamma_steps

        self.cutoff_energy: float = cutoff_energy


    def __repr__(self) -> str:
        return f"{self.gamma_steps} spectra ranging from gamma={self.gamma_min} to gamma={self.gamma_max}."


    def set_spectral_indices(self, gamma_min:float, gamma_max:float) -> None:
        """
        Sets the limits of spectra. gamma_min < gamma_max
        """
        self.gamma_min = gamma_min
        self.gamma_max = gamma_max


    def produce_power_law_spectra(self, response_df=None, energy_grid:np.ndarray=None, cutoff_energy:float=None) -> None:
        """
        Produces a list of spectra, that are needed for bowtie analysis.

        Saves a list of dictionaries containing values for each spectrum to a class attribute "power_law_spectra".

        Parameters:
        -----------
        response_df : {pd.DataFrame} Dataframe containing the response function data.
        energy_grid : {np.ndarray} Array containing the midpoints of the energy grid.
        cutoff_energy : {float} Cutoff energy in MeV. If None, the class attribute "cutoff_energy" is used.
        """

        # The incident energies are needed in a specific format, which is taken care of here.
        validate.validate_response_df_and_grid(response_df=response_df, energy_grid=energy_grid)

        if response_df is not None:
            response_matrix = bowtie_util.assemble_response_matrix(response_df=response_df)
            energy_grid = response_matrix[0]["grid"]["midpt"]
        else:
            energy_grid = energy_grid

        # Check if cutoff energy was provided, if not -> use the class attribute.
        if cutoff_energy is None:
            cutoff_energy = self.cutoff_energy

        # Generates the power law spectra with exponential cutoff, if cutoff_energy > 0, otherwise just power law spectra.
        if cutoff_energy > 0:
            power_law_spectra: list = bowtie_calc.generate_exppowlaw_spectra(energy_grid=energy_grid, 
                                                                             gamma_pow_min=self.gamma_min, gamma_pow_max=self.gamma_max, 
                                                                             num_steps=self.gamma_steps,
                                                                             cutoff_energy=cutoff_energy)
        else:
            power_law_spectra: list = bowtie_calc.generate_pwlaw_spectra(energy_grid=energy_grid,
                                                                         gamma_pow_min=self.gamma_min,
                                                                         gamma_pow_max=self.gamma_max,
                                                                         num_steps=self.gamma_steps)

        # Save the produced power law spectra to class attribute for easy access later.
        self.power_law_spectra: list = power_law_spectra


    def produce_integral_power_law_spectra(self, energy_grid:np.ndarray) -> None:
        """
        Produces a list of spectra, that are needed for bowtie analysis.

        Saves a list of dictionaries containing values for each spectrum to a class attribute "power_law_spectra".
        """

        # Generates the power law spectra
        integral_spectra: list = bowtie_calc.generate_integral_pwlaw_spectra(energy_grid=energy_grid, 
                                                                    gamma_pow_min=self.gamma_min,
                                                                    gamma_pow_max=self.gamma_max, 
                                                                    num_steps=self.gamma_steps)

        # Save the produced power law spectra to class attribute for easy access later.
        self.integral_spectra: list = integral_spectra

