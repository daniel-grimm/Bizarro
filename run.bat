@echo off
REM			the arguments that you can pass are in the following order:
REM			1. number of templates - DEFAULT: 12
REM			2. number of comics (0 - 100 unless you download more) - DEFAULT: 100
REM			3. template scaling factor (finer-grained should probably have more scaling iterations) - DEFAULT: 0.98
REM			4. template scaling iterations - DEFAULT: 7
REM			5. degrees of rotation (the template is rotated 360 / dor times, make sure that 360 % dor == 0) - DEFAULT: 90
REM			6. threshold value (max correlation value - try 0.5 and work from there) - DEFAULT: 0.5
REM			7. max threshold (added to threshold value, indicates when the program should stop trying to match) - DEFAULT: 0.1
REM			8. template penalty (whether the correlation of a match should be penalized depending on the current size of the
REM								template vs the original size of the template) - DEFAULT: ??? (true/false)
REM			9. template penalty modifier (a constant that is added to the ratio depending on the current scaling iteration,
REM											should be a small value like 0.01) - DEFAULT: 0.01
REM			10. debug (inputs: "true" or "false") - DEFAULT: false

SET numberOfTemplates=12
SET numberOfComics=100
SET scalingFactor=0.8
SET scalingIterations=16
SET degreesOfRotation=30
SET thresholdValue=0.6
SET maxThreshold=0.2
SET templatePenalty=true
SET templatePenaltyModifier=0.01
SET debug=false

Final_Project.exe %numberOfTemplates% %numberOfComics% %scalingFactor% %scalingIterations% %degreesOfRotation% %thresholdValue% %maxThreshold% %templatePenalty% %templatePenaltyModifier% %debug%