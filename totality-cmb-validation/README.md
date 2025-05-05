
# Totality Field Theory - CMB High-ℓ Validation

This repository contains the official implementation of the results reported in:

**Section 4 Empirical Validation I: CMB High-ℓ Structure**  
from the Totality Field Theory paper.

It compares the ΛCDM model with the Totality model using Planck 2018 CMB TT power spectrum data, employing harmonic modulations and information-theoretic model comparison (χ², AIC, BIC, p-value).

## 📂 Contents

- `cmb_tt_totality_allN_plot.py`: Main Python script (N=1~5 fitting)
- `data/`: Planck power spectrum data
- `outputs/`: Resulting plot image

## 📚 Data Source

The Planck 2018 CMB TT power spectrum data used in this project is publicly available from the ESA Planck Legacy Archive:

- [Main Archive](https://pla.esac.esa.int/#cosmology)
- [Direct file link]([https://pla.esac.esa.int/pla/slave/API/8.0/cosmology/CMB_spectrum/COM_PowerSpect_CMB-TT-full_R3.01.txt](http://pla.esac.esa.int/pla/aio/product-action?COSMOLOGY.FILE_ID=COM_PowerSpect_CMB-TT-full_R3.01.txt))

File: `COM_PowerSpect_CMB-TT-full_R3.01.txt`


## ▶️ How to Run

1. Install requirements:

```bash
pip install -r requirements.txt
