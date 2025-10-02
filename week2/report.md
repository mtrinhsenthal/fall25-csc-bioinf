# Steps

1. Port matrix functions first - as other files depend on it - and test individual functions
2. Port minimal.py and thresholds.py, and test functions
3. Convert \_pwm.c to Codon
4. Port motif class
5. Write tests for Codon in tests_motif.py - including a custom test runner

# Gotchas

## Issues with subclassing

- We had some issues with subclassing in Codon. The way GenericPositionMatrix is defined in biopython
  GenericPositionMatrix is defined as a subclass of the dictionary class in biopython. Then, the other classes defined in matrix.py were defined as subclasses of the GenericPositionMatrix class. However, even after we defined **getitem** method for the GenericPositionMatrix class, it was not possible for us to index within the subclasses.

  To bypass this, we instead did not define the GenericPositionMatrix class as a subclass of dictionary, but added an attribute to it called data that was a dictionary. We then implemented the **getitem** function for the GenericPositionMatrix such that it would index this data attribute. This way, subclasses of the GenericPositionMatrix were able to call the parent's **getitem** method without any issues.

## Motif.format method

- We were unable to simply import and call the Biopython write methods for the different formats. So instead we had to port them into Codon and put them into the utilities.py file so the program would compile and run correctly.

## parse method (inside init.py)

- We got the error 'File' object has no attribute '**to_py**' when trying to import and call the biopython read methods for the different formats (i.e., without porting them), as instructed, we simply removed all formats that we were not asked to port, so only the read function within minimal.py is called from within this function.

## Tests

- We skip the python test for test_format_clusterbuster (codon test passes); it tries to access Motif.weight which is not possible (problem with biopython?)

- When testing test_pwm_getitem, we were unable to execute some tests that the original biopython tests included, for example, the following code
  ```
    t = counts[i2, i3:i12:i2]
    tc.assertAlmostEqual(t[i0], 0.0)
  ```
  (where i2, i3, i12, i2, and i0 are all integers)
  since t is of type dict[str, list[float]] and codon does not allow to index dictionaries of this type with integers.

# Time Spent
3 - 6 hours every day since Monday Sep 22 (9 days)
