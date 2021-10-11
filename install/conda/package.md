# Create and upload conda package

## Create package

    cd install/conda/
    source ~/miniconda3/bin/activate
    conda build purge
    conda-build -c bioconda -c conda-forge .

## Test package locally

    conda create --name test_package
    conda activate test_package
    conda install -c bioconda -c conda-forge --use-local anacore=2.10.0
    python
    # In interpreter
      from anacore.vcf import VCFIO
      from anacore.db.homo_sapiens.accession import AssemblyAccession
      quit()
    conda deactivate
    conda env remove --name test_package

## Upload to anacoda cloud

    anaconda upload ~/miniconda3/conda-bld/noarch/anacore-*.tar.bz2
    rm ~/miniconda3/conda-bld/noarch/anacore-*.tar.bz2
