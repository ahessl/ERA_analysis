# ERA_analysis

Tools to compare a time series (here tree ring chronology) with climate field data (netCDF), including correlation, compositing, extracting a time series from gridbox.

## Current Tasks
## Current Files
#### *universal_netCDF_script.R - spatial-field correlation of climate data with (tree ring) proxy indices*
##### Functioning Script
  1. Bringing in climate data as a raster
    - crop to user-defined extent
    - temporal crop
    - rename layers for season (three month functioning properly)
       - apply Schulman shift
  2. Associate tree ring data with climate data
    - crop to extent of climate data
  3. Correlations
    - apply correlation to data
    - crop correlations by p-value
    - display map with cropped correlations

##### Debugging Script
  1. detrending tree ring
    - appropriate methods [ discuss additions from AH to understand ]
  2. detrending climate data

##### Test Script
  1. Season/Year-end/Shulman code for universal usage
    - May need to add in flexibility for same year, temporally dislocated proxy months
    - Seasonal conventions/ideas SC not familiar with

#### *Composite.R - seasonal composite of climate variables *
##### Functioning Script
##### Debugging Script
##### Test Script

## Current Overall Tasks

### *Needs*
#### *Needs*
- Plotting of some anomalies has an asymetric range about 0 resulting in plots with colors at "0" (where there should be white). Not sure if we can remedy this. Ideas?

## Extract functions from universalNetCDF to:
### 1. seasNm.R - Seasonalize gridded climate data (netCDF)
- Generic version of plotting similar to the code for universal test #1

## 'Package Components'
functions extracted from code to create generic operations.

### Extract functions from universalNetCDF to:

##### 1. seasNm.R - Seasonalize gridded climate data (netCDF)
  - Four arguments passed to seasNm
    * climDat: climate data
    * SchulmanShift: adjust for seasonal offset in SH. Default is FALSE, calendar year and seasons match
    * lg: option of no lag and one year lag for now.
    * FUN: a function to apply during the compiling of the final product. Accepts regular functions (e.g. mean, sum)
  - created defaults for SchulmanShift and lg for ease of use.
    * SchulmanShift: default FALSE
    * lg: default no lag 0

  #### *Needs*
###### *Needs*
  - Argument defaults  
    * FUN: mean?
  - create more flexibility for different combinations e.g.:
      * Warm vs. Cool season
      * Two-month
      * Specific start and end months --YES!  The ice data are all over the place in terms of season:
      ice accumulation - same year, J-Dec
      seas salts - same year, MJJA-ish
      isotopes - same year, NDJFM-ish (and no Shulman shift)
      MSA - same year, DJFM (no Shulman shift)
      The trick will be to have some clear language!!!


### 2. fullCorr.R - Correlate seasonal climate grids with proxy time series (indices)

  #### *Needs*
  - Holding matrices (Cor, CorT, temp) added into function
  - More generic input of time series to allow for other kinds of data and data formats (ice core data!)

### 3. detrCL.R - Detrend both seasonalized cliamte data and tree ring time series to remove longterm trends

 #### *Needs*
  - default for mod???; should be no detrending
  - Tree rings added
  - This function may not be feasible
  - Detrending is computationally slow - esp ffcsaps (splines)

### 4. compCalc.R - Create composites of climate data based on user-defined quantile of proxy time series
  - Ability to choose between upper ("u"), lower ("l"), or both ("b") when calling the function. No need to assign the function to a variable, it self-produces the results.
  - Produces a composite summary (ComSummary) when answering 'yes' to the prompt in the console. It consists of: #Quantile years, #Mean years, #of each season Quantile, #of each season Mean.

#### *Needs*
  - Error message for uneven seasons (is this necessary now with the ComSummary output?)
  - Let's have a report of the years that were chosen and # for compositing

### 5. ncdfRead.R - imports netcdf file, allows operator to choose variable in the console. Extracts the variable data and determines whether the naming system is "months" or "days and names layers appropriately. Lastly, rotates coordinate system of raster if longitude is not -180 - 180.
* Now searches for lon within dimension title - should fix the issue.
* SC: I think slowness with generating dat needs to be taken as-is for now. Unless we want to dive into paralell processing - which I'm trying to avoid right now....
 
####  *Needs*
  - Tested
  - Clean-up of time units
    * Only have to worry about months and days? or are there other variants? (so far I have only had a problem with "months since...")
  - Runs a little slower than Shawn would like because of rotate, but it may not be avoidable. (Yes, not sure how to avoid, could be temporally cropped first but can't be cropped by spatial extent until units are correct - unless you crop by latitude which is fine)
    * Combine cropping by extent for this portion?
  - Rotate questions
    * Will values always be < 181?: (I notice that all climate data sets I looked at had 0-359 for lon.  It seems that the if/else is not working because of the label calling for longitude is sometimes "lon", sometimes "longitude" etc, I replaced with getvar, but still not working!
    
### 6. ExtractTS
  * A way to do it, some data wrangling behind the scenes will be required and might be very messy to begin with.  
    
#### *Needs*
   - Nothing going yet, idea is to interactively select grid box from plots to extract a mean time series of months that can be sent to treeclim for monthly analysis.
    identifyRaster to get the extent
    gbox <- extent (140, 150, -42, -40)
    gbox_m <- extract(raster, gbox, fun="mean", method="simple")


Build a package called SpatCorr












## Getting Started

Not ready yet

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you have to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [R](https://www.r-project.org/)
* [dplR](https://cran.r-project.org/web/packages/dplR/index.html)



## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Shawn Cockrell** - *Most work* -
* **Amy Hessl** - *Some guidance* -


See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
