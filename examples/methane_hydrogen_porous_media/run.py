"""Calculate diffusion coefficients for methane (CH4)/H2 in porous media"""


def main():
    from porous_media.parameters import FluidMixture
    void_fraction, d_pore, tortuosity = 0.35, 1.e-8, 1.5
    names = ['CH4', 'H2']  # must match what in parameter file!
    molecular_weights = [12.011+1.008*4, 1.008*2]
    from raw_data.read_data import get_LJ_params, read_csv
    LJ_data = read_csv('../../raw_data/LJparams.csv')
    sigma, epsilon = get_LJ_params(names, LJ_data)
    cls = FluidMixture(
        names, molecular_weights, sigma, epsilon, void_fraction, d_pore, tortuosity
    )
    cls.write_params('input.csv')
    cls.write_calculations('output.csv', 300., 101325.)


if __name__ == '__main__':
    main()