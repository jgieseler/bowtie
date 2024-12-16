import filecmp
import math
import matplotlib
import os
import pytest
import bowtie as bow
import numpy as np
import pandas as pd


"""
Install dependencies for tests:
pip install flake8 pytest pytest-doctestplus pytest-cov pytest-mpl

To create/update the baseline images, run the following command from the base package dir:
pytest --mpl-generate-path=tests/baseline bowtie/tests/test.py

To run the tests locally, go to the base directory of the repository and run:
pytest -rP --mpl --mpl-baseline-path=baseline --mpl-baseline-relative --mpl-generate-summary=html --cov=bowtie bowtie/tests/test.py
"""


@pytest.mark.mpl_image_compare(remove_text=True, deterministic=True)
def test_bowtie():
    response_df = pd.read_csv("sixs_side0_electron_responses.csv", index_col="incident_energy")
    #
    energy_min = 0.01
    energy_max = 50.
    bowtie = bow.Bowtie(energy_min=energy_min, energy_max=energy_max, data=response_df)
    #
    gamma_min = -5.5
    gamma_max = -2.5
    num_of_spectra = 100
    spectra = bow.Spectra(gamma_min=gamma_min, gamma_max=gamma_max, gamma_steps=num_of_spectra)
    #
    spectra.produce_power_law_spectra(response_df)
    #
    channel = "E1"
    e1_results = bowtie.bowtie_analysis(channel=channel, spectra=spectra, plot=False)

    assert e1_results['geometric_factor'] == 0.00046108443865681666
    assert e1_results['geometric_factor_errors'] == {'gfup': np.float64(1.211669671885579e-05), 'gflo': np.float64(8.045747419734267e-06)}
    assert e1_results['effective_energy'] == np.float64(0.0607020174186593)
    # assert e1_results['fig'] == 
    # assert e1_results['axes'] == 

    new_gamma_min = -3.5
    new_gamma_max = -1.5
    spectra.set_spectral_indices(gamma_min=new_gamma_min, gamma_max=new_gamma_max)
    spectra.produce_power_law_spectra(response_df=response_df)

    new_e1_results = bowtie.bowtie_analysis(channel=channel, spectra=spectra, plot=False)

    assert new_e1_results['geometric_factor'] == 0.000744084026271278

    assert math.isclose(new_e1_results['geometric_factor_errors']['gfup'], np.float64(4.407576153642298e-05))
    assert math.isclose(new_e1_results['geometric_factor_errors']['gflo'], np.float64(3.07655324791197e-05))
    assert math.isclose(new_e1_results['effective_energy'], np.float64(0.0713699882077903))
    # assert new_e1_results['fig'] == 
    # assert new_e1_results['axes'] ==

    all_channels_results = bowtie.bowtie_analysis_full_stack(spectra=spectra, plot=True)

    assert len(all_channels_results) == response_df.shape[1]

    # only check the last result
    assert all_channels_results[6]['geometric_factor'] == 0.2955529518148398
    # assert all_channels_results[6]['geometric_factor_errors'] == {'gfup': np.float64(0.027996508853041446), 'gflo': np.float64(0.01830299979208344)}
    assert math.isclose(all_channels_results[6]['geometric_factor_errors']['gfup'], np.float64(0.027996508853041446))
    assert math.isclose(all_channels_results[6]['geometric_factor_errors']['gflo'], np.float64(0.01830299979208344))
    assert math.isclose(all_channels_results[6]['effective_energy'], np.float64(7.0097596149373445))
    # assert all_channels_results[6]['fig'] == 
    # assert all_channels_results[6]['axes'] ==

    # test saving
    filename = "test.csv"
    column_names = response_df.columns
    bow.bowtie_util.save_results(results=all_channels_results, filename=filename, column_names=column_names, save_figures=True)

    filecmp.cmp('test.csv', 'bowtie/tests/test_org.csv')

    for i in range(1,8):
        assert os.path.exists(f'E{i}_bowtie.png')

    # return last produced fig for mpl pytest
    return all_channels_results[6]['fig']
