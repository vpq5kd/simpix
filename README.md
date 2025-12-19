# simpix

C++ starter code
* simpix_start.cpp
use make to build this example

Usage: simapix_start image1 image2 <output=out.png>

Python starter code
* simpix_start.py

Usage: simapix_start image1 image2 <output=out.png>

***SIMPIX SOLUTIONS AND WRITEUP***
Program is unchanged from the initial starter code and can be run and built the exact same way. That being said, you are able to change the initial temperature, melting iterations, and number of iterations per temp step in the code itself. I, however, would not recommned doing this as the current parameters appear to be working *quite* well. 

**Note:** Unlike my salesman solution, this solution is still written with linear temperature decreasing, static internal iterations, and a fixed start temperature. Though this algorithm could be changed, the current code outputs a desirable solution in an extremely low amount of time, so it is uneccessary to change it at this time.

**640x180** "Frisbe to Rotunda" occured in 8.96s. Its out-file is located in rotundaOutFile.png and its progression collage is located in scottToRotunda.png

**1920x1080** "Pollock to Pollock" occured in 10.56s. Its out-file is located in pollockOut.png and its progression collage is located in pollockToPollock.png
