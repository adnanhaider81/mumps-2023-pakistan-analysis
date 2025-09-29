# Contributing

Thanks for your interest in improving this analysis. Follow the steps below.

## Ground rules
- Protect privacy. Do not commit patient data or any restricted dataset.
- Keep methods reproducible. Any change in parameters should be in `config/config.yaml`.
- Match tool versions to the `env/environment.yml` when possible.

## Workflow
1. Fork the repo then create a feature branch.
2. Make focused changes with clear commit messages.
3. Run a dry run before opening a PR:
   ```bash
   conda activate mumps-2023-env
   export NCBI_EMAIL="your@email"
   snakemake -s workflow/Snakefile -n -c 4
   ```
4. Verify `make test` works.
5. Open a pull request that describes the change, the motivation, and any parameter updates.

## Code style
- Python scripts live in `analysis/scripts`. Keep functions small. Prefer argparse. Write helpful `--help` text.
- Use TSV for tables and FASTA for sequences. Put large outputs under `results/`.
- Avoid hard coded paths. Use command line flags and config.

## Data and secrets
- Do not check in raw FASTQ or BAM files. Use `data-private/` locally.
- Set `NCBI_EMAIL` and optional `NCBI_API_KEY` in your shell only.

## Questions
Open an issue with a short description and steps to reproduce if relevant.
