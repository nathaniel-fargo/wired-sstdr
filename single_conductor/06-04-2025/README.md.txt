# Single Conductor Tests

Today we finished the setup for the single conductor experiments.

We're taking measurements at different distances away.

We'll make a CSV file that records the IDs for the measurements and their corresponding parameters. 

The single conductor was 40' and we took off 33" which makes it 37'3"

The setup is WILMA -> D00 (BNC cable) -> BNC to alligator clips -> Single conductor over ground plane (37'3") at variable distance

The Z00 measurement was from the WILMA -> D00 (BNC Cable) -> BNC to alligator clips

All of the other test data is labeled T[00-99]_[0-9]. 
For now, the first number seems to line up with the number of foam blocks placed in between, so the depth can be approximated as (1.5 + 1.1 * n), with an initial depth of 1.5 and the foam blocks are a little more than an inch tall. 
The second number indicates which test iteration it was. 