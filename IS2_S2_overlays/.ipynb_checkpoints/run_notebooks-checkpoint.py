import papermill as pm

example = 'P1'
beam_name = 'gt1l'

# Define the list of parameters
params_list = [
    {"example": "P1", "beam_name": "gt1l"},
    {"example": "P1", "beam_name": "gt2l"},
    {"example": "P1", "beam_name": "gt3l"},
    {"example": "P1", "beam_name": "gt1r"},
    {"example": "P1", "beam_name": "gt2r"},
    {"example": "P1", "beam_name": "gt3r"},
]

# Iterate over the parameters and execute the notebook
for i, params in enumerate(params_list):
    pm.execute_notebook(
        '/home/jovyan/GitHub/ICESat-2-sea-ice-tools/IS2_S2_overlays/plot_s2_atl07_cloud.ipynb',
        f'/home/jovyan/GitHub/ICESat-2-sea-ice-tools/IS2_S2_overlays/output_notebook_{i}.ipynb',
        parameters=params
    )
