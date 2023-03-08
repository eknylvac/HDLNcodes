# HDLNcodes
In this repository, you can find the MATLAB codes of 'Design of Urban Logistics Networks to Support Home Delivery' paper. Besides the functions that were explained below, you can find the main code of our design and other minor functions that were used in this design.

## Contents
* [DemandWeight](https://github.com/eknylvac/HDLNcodes/blob/main/demandweight.m)
* [LocateDC](https://github.com/eknylvac/HDLNcodes/blob/main/locateDC.m)
* [LinehaulLoadDemand](https://github.com/eknylvac/HDLNcodes/blob/main/linehaulloaddemand.m)

### DemandWeight
This function utilizes the proximtiy factor algorithm to find the percentage of interactions between nodes in space using their weight and distance data. In logistics network design, this weight is generally the population data. Inputs are as follows:
* w: weight.
* D: distance.
* in: input structure. We use the inbound distance (dinb) to a DC in this input structure.

Outputs are as follows:
* W: the matrix that shows the percentages of interactions between nodes. The sum of all values in this matrix adds up to one.
* pxf: calculated proximitf factor using the inbound distance. This value is the optimal proximity factor for our design. MATLAB's fminsearch optimization algorithm is used for the optimization.

For more information on the proximity factor please refer to
Yalvac, E., & Kay, M. G. (2022). [Synthetic demand flow generation using the proximity factor](https://assets.researchsquare.com/files/rs-1918195/v1_covered.pdf?c=1659967372).


### LocateDC
This function finds the locations of distribution centers (DCs) by allocating a point's demand to the nearest DC until it reaches its capacity. Inputs are as follows:
* f: demand vector of aggregate demand points (ADPs).
* P: coordinates of ADPs.
* a: area vector of ADPs (to use it in the aggregate distance calculation).
* in: input structure. We use, fmax, maximum load per hour at a DC.

Outputs are as follows:
* X: coordinate matrix of DCs
* F: a matrix that shows how much a node is served by which DC.
* TC: Latest iteration's total cost.
* TC0: Lowest total cost of all iterations.

For more information on the aggregate distance please refer to
Yalvac, E., & Kay, M. G. (2022). [Synthetic demand flow generation using the proximity factor](https://assets.researchsquare.com/files/rs-1918195/v1_covered.pdf?c=1659967372).

### LinehaulLoadDemand
This function finds the load demand of linehaul for each DC. Inputs are as follows:
* X: coordinates of DCs.
* W: the matrix that shows the percentages of interactions between DCs.
* Q: total population.
* a: area vector of DCs (to use it in the aggregate distance calculation).
* in: input structure. We are using multiple parameters that were stored in this structure.

Outputs are as follows:
* fH: linehaul between DCs in terms of load per hour.
* T: travel plus loading and unloading time.
* Aw: the matrix that shows the percentages of interactions between DCs except the diagonal is zero.

For more information on the aggregate distance please refer to
Yalvac, E., & Kay, M. G. (2022). [Synthetic demand flow generation using the proximity factor](https://assets.researchsquare.com/files/rs-1918195/v1_covered.pdf?c=1659967372).
