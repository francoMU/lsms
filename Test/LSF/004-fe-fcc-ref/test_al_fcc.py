import contextlib
import os
import shutil
import unittest
from pathlib import Path

from gfalm import main
from gfalm.utils import read_input_file
from gfalm.wrap import output_data_mod


@contextlib.contextmanager
def change_dir():
    prev_dir = os.getcwd()
    try:
        os.chdir(Path(__file__).parent)
        yield
    except Exception as e:
        os.chdir(prev_dir)
        raise e


class TestAlFcc(unittest.TestCase):

    def setUp(self) -> None:

        try:
            shutil.rmtree('struct', ignore_errors=True)
            shutil.rmtree('job.test-fcc', ignore_errors=True)
        except FileNotFoundError:
            pass

    def test_run_al_fcc(self):
        """
        Test the output of one big iteration for fcc-Al.
        """

        with change_dir():
            input_data = read_input_file('al_fcc.in')
            #
            # This is the first test of the series: remove the cached SCMs
            # if they have already been created
            #

            #
            # Run the main code
            #
            return_state = main.run(**input_data)

            #
            # Return status 1 tells that the total energy did not converge
            #
            self.assertEqual(return_state.stat, 1)

            latt_site = output_data_mod.get_lattice_site_output(1)

            bgfm_energy = -483.969255

            #
            # Renormalized values
            #
            bgfm_Qlm1 = 3.0  # l = 0, m = 0
            bgfm_Qlm2 = -0.0028263464304121  # l = 4, m = 0
            bgfm_Qlm3 = -0.0023886987109995  # l = 4, m = 4

            self.assertAlmostEqual(latt_site.energy, bgfm_energy, delta=1e-6)
            self.assertAlmostEqual(latt_site.qlm_av[0], bgfm_Qlm1, delta=1e-8)
            self.assertAlmostEqual(latt_site.qlm_av[20], bgfm_Qlm2, delta=1e-8)
            self.assertAlmostEqual(latt_site.qlm_av[24], bgfm_Qlm3, delta=1e-8)


if __name__ == '__main__':
    unittest.main(verbosity=2, buffer=False)
