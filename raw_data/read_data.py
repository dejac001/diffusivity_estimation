def read_csv(f_name):
    with open(f_name) as f:
        header = list(map(eval, next(f).rstrip('\n').split(',')))
        data = {key: [] for key in header}
        for line in f:
            for key, val in zip(header, line.rstrip('\n').split(',')):
                try:
                    data[key].append(float(val))
                except ValueError:
                    data[key].append(val)
    return data


def get_LJ_params(components, LJ_data):
    """Get LJ params from csv file

    :param components:
    :return: LJ params sigma and epsilon
    """
    sigma = []
    epsilon = []
    for key in components:
        assert key in LJ_data['Substance'], 'No parameters found for {} in LJparams.csv'.format(key)
        index = LJ_data['Substance'].index(key)
        sigma.append(LJ_data['sigma [angstrom]'][index])
        epsilon.append(LJ_data['epsilon [K]'][index])
    return sigma, epsilon
