# ERA_analysis

Tools to compare a time series (here tree ring chronology) with climate field data (netCDF), including correlation, compositing, extracting a time series from gridbox.

## Current Tasks

*Needs*
Plotting of some anomalies has an asymetric range about 0 resulting in plots with colors at "0" (where there should be white). Not sure if we can remedy this. Ideas?

Extract functions from universalNetCDF to:
1. seasNm.R - Seasonalize gridded climate data (netCDF)

  *Needs*
  - create defaults for hem, lg, FUN for quick processing
    hem: default NH
    lg: default no lag 0
    FUN: default....mean?
  - create more flexibility for different combinations e.g.:
      * Warm vs. Cool season
      * Two-month
      * Specific start and end months --YES!  The ice data are all over the place in terms of season:
      ice accumulation - same year, J-Dec
      seas salts - same year, MJJA-ish
      isotopes - same year, NDJFM-ish (and no Shulman shift)
      MSA - same year, DJFM (no Shulman shift)
      The trick will be to have some clear language!!!


2. fullCorr.R - Correlate seasonal climate grids with proxy time series (indices)

  *Needs*
  - Holding matrices (Cor, CorT, temp) added into function
  - More generic input of time series to allow for other kinds of data and data formats (ice core data!)

3. detrCL.R - Detrend both seasonalized cliamte data and tree ring time series to remove longterm trends

  *Needs*
  - default for mod; should be no detrending
  - Tree rings added
  - This function may not be feasible

4. compCalc.R - Create composites of climate data based on user-defined quantile of proxy time series

  *Needs*
  - Error message for uneven seasons
  - What is the extra argument on line 115 for com.l ("5")???
  - Let's have a report of the years that were chosen and # for compositing

5. ncdfRead.R - imports netcdf file, allows operator to choose variable in the console. Extracts the variable data and determines whether the naming system is "months" or "days and names layers appropriately. Lastly, rotates coordinate system of raster if longitude is not -180 - 180.
 
  *Needs*
  - Tested
  - Clean-up of time units
    * Only have to worry about months and days? or are there other variants? (so far I have only had a problem with "months since...")
  - Runs a little slower than Shawn would like because of rotate, but it may not be avoidable. (Yes, not sure how to avoid, could be temporally cropped first but can't be cropped by spatial extent until units are correct - unless you crop by latitude which is fine)
    * Combine cropping by extent for this portion?
  - Rotate questions
    * Will values always be < 181?: (I notice that all climate data sets I looked at had 0-359 for lon.  It seems that the if/else is not working because of the label calling for longitude is sometimes "lon", sometimes "longitude" etc, I replaced with getvar, but still not working!
    
5. ExtractTS
    *Needs*
    --Nothing going yet, idea is to interactively select grid box from plots to extract a mean time series.
    identifyRaster
    extract (raster, coords)

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
