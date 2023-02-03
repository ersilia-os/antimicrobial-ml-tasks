#!/bin/bash
# Call the zairachem command, redirect logs, and write a line into "runs.csv"
#
# Requirements:
#     - ZairaChem must be installed - therefore conda environment "zairachem" exists.
#     - Subdirectories "input", "model", "test" and "log" must exist in the current directory.
#     - Depending on the command, the files "input/input.csv", "input/train.csv" or "input/test.csv" must exist.
#
# Examples of use:
#     bash test.sh split similarity
#     bash test.sh fit
#     bash test.sh fit --lazy
#     bash test.sh predict

command="$1"
option="$2"

eval DIR_MODEL_RUNS_LOG="~/models"  # Where to store the runs log
VALID_COMMANDS='["split", "fit", "predict"]'

if [[ $VALID_COMMANDS =~ "\"$command\"" ]]
then
    # Command is valid
    printf -v ts_start '%(%Y%m%d_%H%M)T' -1  # Get current datetime in yymmdd_hhmm format
    echo "Calling zairachem $command - $ts_start"
    pwd

    if [ "$command" = "fit" ]
    then
        conda run -n zairachem zairachem fit -i input/train.csv -m model $option --clean >log/fit_$ts_start.stdout.txt 2>log/fit_$ts_start.log
        exit_code=$?
    elif [ "$command" = "predict" ]
    then
        conda run -n zairachem zairachem predict -i input/test.csv -m model -o test >log/predict_$ts_start.stdout.txt 2>log/predict_$ts_start.log
        exit_code=$?
    elif [ "$command" = "split" ]
    then
        conda run -n zairachem zairachem split -i input/input.csv --split_criterion $option >log/split_$ts_start.stdout.txt 2>log/split_$ts_start.log
        exit_code=$?
    fi

    printf -v ts_end '%(%Y%m%d_%H%M)T' -1  # Get current datetime in yymmdd_hhmm format
    echo "Finished zairachem $command - exit code $exit_code - $ts_end"
    # Write a line in runs log
    mkdir -p $DIR_MODEL_RUNS_LOG
    echo `pwd`,$command,$ts_start,$ts_end,$exit_code >>$DIR_MODEL_RUNS_LOG/runs.csv
else
    # Command is not valid
    echo "Error - the valid commands are $VALID_COMMANDS"
fi
