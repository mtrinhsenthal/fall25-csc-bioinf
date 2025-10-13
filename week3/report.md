
## Steps
- copy and paste the relevant files for these 3 tests from the biotite repo
- convert from cython to codon by removind cdefs, rearranging type annotations, and adding type annotations
- try running codon, go through the errors, and correct one by one (mostly about fixing type annotations)
- make duplicate files for python version (since python complains about import statements in codon)
- add timing functions to both, create evaluate.sh that calls both test files

## Gotchas
The codon tests seem to take longer than the python tests, though I am not sure why.

## Timem Estimate
Time estimate for this assignment: 7 hours