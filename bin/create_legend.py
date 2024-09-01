import matplotlib.pyplot as plt

def generate_legend(symbol_description_fillcolor_pairs, output_file='legend-0.png'):
    fig, ax = plt.subplots()

    for i, (symbol, description, fill, color, size) in enumerate(symbol_description_fillcolor_pairs):
        marker_params = {'marker': symbol, 'label': description, 'color': color, 's': size}

        if not fill:
            marker_params['facecolor'] = 'none'  # empty marker

        ax.scatter([], [], **marker_params)

    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Legend')

    # Remove x and y ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Remove x and y axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.5, dpi=300)

symbol_description_fillcolor_pairs = [
    ((12,1,0), 'Blazar candidate - HBL type', True, '#ff9933',120),
    ((12,1,0), 'Blazar candidate - IBL type', True, '#66ccff',120),
    ((12,1,0), 'Blazar candidate - LBL type', True, '#0033cc',120),
    ((12,1,0), 'Radio emitting Seyfert/QSO candidate ', True, '#555566',120),
    ('o', 'X-ray source possible blazar candidate', False, 'blue',60),
    ('o', 'Radio source possible blazar candidate', True, '#ff0066',60),
    ((8,1,0), 'Radio quiet AGN candidate ', True, '#ff0066',100),
    ((4,1,0), 'Star X-ray counterpart candidate', True, '#00cc66',100),
    ('$Clu$', 'Cluster of galaxies', True, '#993300',180),
    ('o', 'Radio source', False, '#ff0066',60),
    ((4,1,0), 'Optical source', True, 'gray',80),
    ('x', 'X-ray source', True, 'blue',60),
    ('^', 'Gamma-ray source', False, '#7e0de8',100),
    ('o', 'Fermi 4LAC source', False, '#ff1aff',140),
    ('o', 'IHBL source', False, 'green',120),
    ('o', 'QSO from milliquas catalog', False, '#84e184',60),
    ('*', '3HSP source', False, '#ffcc99',120),
    ('D', '5BZcat source', True, '#ffcc99',40),
    ('s', 'CRATES source', False, '#6267ca',80),
    ('p', 'Pulsar', True, '#8000ff',100),
    ('H', 'GRB', False, '#555555',100),
    ('$CV$', 'Cataclysmic variable', True, '#ff0066',120),
    ('$MC$', 'Molecular Cloud', True, '#ff0066',150),
    ('$SC$', 'Stellar  Cluster', True, '#ff0066',120),
    ('$XRB$', 'X-ray binary', True, '#ff0066',240),
]
generate_legend(symbol_description_fillcolor_pairs)

