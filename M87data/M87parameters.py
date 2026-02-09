M87parameters = {
    "EHT2018": {
        "D_Mpc": 16.8,
        "MBH_MSUN": 6.5e9,
        "theta_view": 17,
        "redshift": 0.00428
    },
    "Lucchini19": {
        "D_Mpc": 16.8,
        "MBH_MSUN": 6.5e9,
        "theta_view": 14,
        "redshift": 0.00428
    }
}

def get_M87_parameters(source="EHT2018"):
    if source in M87parameters:
        return M87parameters[source]["D_Mpc"], \
            M87parameters[source]["MBH_MSUN"], \
            M87parameters[source]["theta_view"], \
            M87parameters[source]["redshift"]
    else:
        print(f"parameters from {source} not known!")
