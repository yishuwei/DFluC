# Program DFluC

Program DFluC is a MATLAB program for Detrended Fluctuation Analysis (DFA).

As compared to other DFA software, Program DFluC is implemented with special attention to 
edge cases such as partly unobserved and/or irregularly sampled data.

Prolegomenon
------
DFA is a method of estimating the Hurst exponent, which is used to quantify self-similarity 
and long-range dependence of a 1-dimensional process.

[Note. For convenience this documentation does not distinguish between a proper process and 
an array of (observed and possibly some unobserved or "missing") values sampled from that 
process. The process/array is denoted by `y` below and also referred to as "data". ("Signal" 
is another common term but is not used in this documentation.)]

DFA consists of the following steps:

1. Divide data into non-overlapping segments or "boxes" of a certain length.
2. Subtract local polynomial trends fit within each box and calculate the root mean square 
of the residuals (root-mean-square fluctuation).
3. Repeat the previous steps for different box sizes.
4. Estimate the power-law relationship between root-mean-square fluctuation and box size 
by means of regression.

When dealing with partly unobserved or irregularly sampled data, Program DFluC fits local 
polynomial trends to data points *in situ* and avoids "stitching" or "warping" data, in order 
not to create artificial jumps or alter the autocorrelation structures. You can of course 
"stitch" and "warp" data yourself, if you want to, by removing all unobserved entries and 
running Program DFluC without providing the actual sample points. You would then need to 
interpret the result with caution.

User Guide
------
You should be able to run Program DFluC in any reasonably recent version of MATLAB. To 
estimate the Hurst exponent `H` of a regularly sampled (i.e., sampled at constant intervals) 
random walk-like process `y`, you can simply execute
```matlab
H = dfa(y)
```
in MATLAB. Program DFluC allows `NaN` entries in `y` and interprets them as unobserved 
(missing) values.

**Important Note.** Given a noise-like process, it is common to estimate the Hurst exponent 
of its (random walk-like) integrated process. You can do so with Program DFluC by executing
```matlab
H = dfa(cumsum(y - mean(y)))
```
in cases where `y` is noise-like and regularly sampled without missing values. The cumulative 
sum `cumsum(y - mean(y))` is often called the "profile" in the literature. Many other DFA 
implementations automatically apply cumulative sum to any data whatsoever. However, this trick 
only makes sense for regularly sampled data without missing values. If this is not the case, 
then you need more clever methods to derive the integrated process (a topic that is beyond the 
scope of this documentation). For regularly sampled data with missing values, a simple 
workaround is to concatenate the observed parts and then take the cumulative sum of the 
resulting "stitched" process. The effect of such a procedure on DFA depends on many factors 
and has been documented in detail by others [Phys. Rev. E 81: 031101]. You should critically 
evaluate whether it is appropriate for your data. 

Program DFluC can be run with up to 5 additional input arguments:
```matlab
H = dfa(y, x, detrend_order, [xstart xend], box_sizes, plotting)
```
(setting any of the additional arguments to `[]` means using the default for that argument), 
where

- `x` is an array indicating the sample points (where the values of `y` are sampled). This is 
helpful if `y` is irregularly sampled, or if you want to specify `[xstart xend]`. If `x` is not 
provided, then `y` is assumed to be regularly sampled.

- `detrend_order` specifies the degree of polynomials to be used for local detrending. Default 
is 2.

- `[xstart xend]` defines the portion of data for which the Hurst exponent is to be estimated. 
`xstart` and `xend` are assumed to have the same unit as `x`. If `[xstart xend]` is not 
provided, then the Hurst exponent will be estimated for the entire range of `x`.

- `box_sizes` is an array indicating several box sizes to be used for local polynomial 
detrending. Each entry is assumed to have the same unit as `x`. If `box_sizes` is not provided, 
then 50 automatically determined box sizes will be used.

  You should choose box sizes such that each box covers a sufficient number of data points, 
while at least a few boxes are available (i.e., not too small, not too big).

- `plotting` is `true` if a plot of log(root-mean-square fluctuation) versus log(box size) is 
to be shown, `false` otherwise. The slope of the least-squares fit line is the estimated Hurst 
exponent `H`. Default is `false`.

  If you could see an apparent crossover in this plot, then it's better to re-estimate the 
Hurst exponent for small and large box sizes separately. Crossovers may be due to finite data 
length, noise/artifact contamination, or intrinsic properties of the process (a complete survey 
of the reasons for crossovers is beyond the scope of this documentation) and thus warrant 
careful study (e.g., by checking whether crossovers depend on `detrend_order`).

The Hurst exponent is a dimensionless measure and does not depend on the units of `x` and `y` 
as long as you make sure that `x`, `[xstart xend]` and `box_sizes` are expressed in the same 
unit if you set them manually. To reduce rounding errors, it is best if you choose a unit 
that renders `x`, `[xstart xend]` and `box_sizes` integer-valued.

License, Disclaimer, and Reference
------
Program DFluC is intended to be an academic software program. Permission to use, copy, modify, 
and distribute the software and its documentation for not-for-profit purposes is granted to 
any person obtaining a copy of the source code, provided that this permission notice appear in 
all copies. For other uses, please contact the author (Y. Wei).

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS 
SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL 
THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY 
DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF 
CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR 
PERFORMANCE OF THIS SOFTWARE.

If you publish results obtained with Program DFluC, you may include the following citation:

Colombo, Wei, Ramautar, Linkenkaer-Hansen, Tagliazucchi, Van Someren. More severe insomnia 
complaints in people with stronger long-range temporal correlations in wake resting-state EEG. 
Frontiers in Physiology 7 (2016). DOI: 10.3389/fphys.2016.00576
