# ERA_analysis

Tools to compare a time series (here tree ring chronology) with climate field data (netCDF), including correlation and compositing.

## Current Tasks

Extract functions from universalNetCDF to:
1. seasNm.R - Seasonalize gridded climate data (netCDF)

*Needs x, y, z*

2. fullCorr.R - Correlate seasonal climate grids with tree ring time series (indices)

*Needs x, y, z*

3. detrCL.R - Detrend both seasonalized cliamte data and tree ring time series to remove longterm trends

*Needs x, y, z*


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
