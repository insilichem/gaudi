case ${CI_OS} in
    windows*)
        eval "$(${CONDA}/condabin/conda.bat shell.bash hook)";;
    macOS*)
        eval "$(${CONDA}/condabin/conda shell.bash hook)";;
    *)
        eval "$(conda shell.bash hook)";;
esac
