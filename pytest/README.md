# Test NanoAODPlus by histograming the extracted data.

The histograms and their expected mean value are defined in a 
yaml file. During the test creation the data file is read and 
histograms are filled.

The implementation follows the example at  https://docs.pytest.org/en/6.2.x/example/nonpython.html

## Test definitions

The tests are defined by a yaml file. Following syntax is supported

    ---
    events: 100
    histo1D:
        Jet_pt:                         # name of histogram and of the variable to be histogramed
            title: "Jet pt"             # title of histogram
            bins: [100, 0, 1000.]       # histogram binning
            test_mean: [30,40]          # resulting histogram has to have a mean in given range

## Installation:
   
Modern ROOT, python3, as provided by docker container rootproject/root:6.24.06-centos7
In addition pytest and pyyaml is required 

## Usage:

    pytest pytest/test_nanoaod.yaml --input=$CMSSW_VERSION/src/test.root --output=nanoaod_histos.root --junit-xml=test_nanoaod.junit.xml

Author: Dietrich Liko \<Dietrich.Liko@oeaw.ac.at\>
