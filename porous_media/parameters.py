import math


class PureFluid:
    """Only one component is moving through porous media"""
    def __init__(self, void_fraction, d_pore, tortuosity, molecular_weight=None):
        """

        :param molecular_weight:            molecular weight [kg/mol] for each component
        :param void_fraction: void fraction of porous media
        :param d_pore:        nominal pore diameter [meters]
        :param tortuosity:    tortuosity of porous media
        """
        self.void_fraction = void_fraction
        self.d_pore = d_pore
        self.molecular_weight = molecular_weight
        self.tortuosity = tortuosity
        self.R = 8.314  # J/mol/K [=] m**3*Pa/mol/K

    def knudsen(self, T):
        """Knudsen diffusivity of component i; derived from kinetic theory of gases

        Units:
            R = m**3*Pa/mol/K [=] m*m*kg/mol/s/s/K

            sqrt(RT/MW) [=] sqrt(m*m/s/s)


        :param i: component name
        :return:  Knudsen diffusivity m^2/s
        """
        return self.void_fraction*self.d_pore/self.tortuosity/3.*math.sqrt(8*self.R*T/math.pi/self.molecular_weight)

    def write_params(self, file_name):
        with open(file_name, 'w') as f:
            for key, val in self.__dict__.items():
                if isinstance(val, float):
                    f.write('{},{},{}\n'.format(key, ' ', val))
                elif isinstance(val, dict):
                    f.write('{},{},{}\n'.format(key, ' ', ' '))
                    for key2, val2 in val.items():
                        if isinstance(key2, tuple):
                            key2 = ' '.join(key2)
                        f.write('{},{},{}\n'.format(' ', key2, val2))
                elif isinstance(val, list):
                    f.write('{},{},{}\n'.format(key, ' ', ' '))
                    for val2 in val:
                        f.write('{},{},{}\n'.format(' ', val2, ' '))
                elif isinstance(val, str) or val is None:
                    f.write('{},{},{}\n'.format(key, ' ', val))
                else:
                    raise Exception('Val type not found for {}:{}'.format(key, val))


class FluidMixture(PureFluid):
    """Multiple components flowing through porous media"""
    def __init__(self, components, molecular_weight, sigma, epsilon_molecular, *args):
        """

        :param components:
        :type components: list
        :param molecular_weight:
        :type molecular_weight: list
        :param sigma:   collision diameter in Angstrom, get from Poling et al "The Properties of gases and liquids"
        :type sigma: list
        :param epsilon_molecular: molecular epsilon in K,get from Poling et al "The Properties of gases and liquids"
        :type epsilon_molecular: list
        :param args: args to pass to super
        """
        PureFluid.__init__(self, *args)
        self.molecular_weight_i = {
            key: val for key, val in zip(components, molecular_weight)
        }
        self.sigma_i = {
            key: val for key, val in zip(components, sigma)
        }
        self.epsilon_i = {
            key: val for key, val in zip(components, epsilon_molecular)
        }
        self.component_i = components

    def knudsen_i(self, i, j, temperature, pressure):
        """

        :param i: component name
        :return: Knudsen diffusivity for component *i*
        """
        self.molecular_weight = self.molecular_weight_i[i]
        return self.knudsen(temperature)

    def sigma_ij_rule(self, i, j):
        """Lennard-Jones combining rule for sigma"""
        return (self.sigma_i[i] + self.sigma_i[j]) / 2.

    def epsilon_ij_rule(self, i, j):
        """Lennard-Jones combining rule for epsilon"""
        return math.sqrt(self.epsilon_i[i] * self.epsilon_i[j])

    def omega_ij(self, w):
        """Bird 1960 p. 746

        :return: Temperature-dependent collision integral (dimensionless)
        """

        A1, B1, A2, B2 = 1.06036, 0.15610, 0.19300, 0.47635
        A3, B3, A4, B4 = 1.03587, 1.52996, 1.76474, 3.89411
        value = (
                A1 / pow(w, B1)
                + A2 / math.exp(B2 * w)
                + A3 / math.exp(B3 * w)
                + A4 / math.exp(B4 * w)
        )
        assert 0.5 < value < 2.7, 'Value predicted of %2.3f is not reasonable' % value
        return value

    def molecular_ij(self, i, j, T, P):
        """Molecular diffusivites estimated by Chapman-Enskog

        :param i: sorbate i name
        :param j: sorbate j name
        :param T: temperature in K
        :param P: pressure in atm
        :return: binary diffusion coefficient (m^2/s)
        """

        # convert molecular weights to g/mol to apply formula
        M_i = self.molecular_weight_i[i]*1000.
        M_j = self.molecular_weight_i[j]*1000.

        return 1.858e-3 * math.sqrt(T*T*T * (1. / M_i + 1. / M_j)) / (
                P * self.sigma_ij_rule(i, j) * self.sigma_ij_rule(i, j) *
                self.omega_ij(T / self.epsilon_ij_rule(i, j))
        )/100./100.

    def effective_macropore_i(self, i, j, temperature, pressure):
        """

        :param temperature: temperature in K
        :param pressure: pressure in Pa
        :param i:       component i name
        :param j:       other_component name
        :return:        effective macropore diffusivity m^2/s
        """
        P_atm = pressure/101325.
        return self.void_fraction/self.tortuosity/(
                1. / self.molecular_ij(i, j, temperature, P_atm) + 1. / self.knudsen_i(i, j, temperature, pressure)
        )

    def write_calculations(self, output_file, temperature, pressure):
        with open(output_file, 'w') as f:
            for i in self.component_i:
                for j in self.component_i:
                    if i == j:
                        continue
                    for attr in [
                        'effective_macropore_i',
                        'knudsen_i',
                        'molecular_ij'
                    ]:
                        func = getattr(self, attr)
                        f.write('%s,%s,%s [m^2/s],%e\n' % (i, j, attr, func(i, j, temperature, pressure)))
