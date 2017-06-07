#!/bin/bash
# VERSION: 1.1
# AUTHORS: Patrice Dehais - SIGENAE Toulouse (sigenae-support@listes.inra.fr)
# COPYRIGHT: 2015 INRA
# LICENSE: GNU GPLv3

function error ()
{
    echo 1>&2 ERROR: $*
    echo 1>&2 "Usage: $0 [qsub options] shell_command_file"
    exit 1
}

function father_call ()
{
    if [ $# -eq 0 ] ; then
	error "missing parameter"
    fi
    cmd_file=`echo $*|awk '{print $NF;}'`
    if [ ! -e $cmd_file -o -d $cmd_file ] ; then
	error "last parameter $cmd_file should be a file"
    fi
    m=`cat $cmd_file|wc -l`
    if [ $? -ne 0 ] ; then
	error "can't retrieve number of lines of file $cmd_file"
    fi
    if [ $m -eq 0 ] ; then
	error "$cmd_file is empty"
    fi

    # store parameters, exept last one
    i=1
    for v in "$@" ; do
        if [ $i -eq $# ] ; then break ; fi
	    ARGV[$i]=$v
	    i=$(($i+1))
    done
    # submit job array
    qsub -t 1-$m -v QARRAY_CMD=$cmd_file "${ARGV[@]}" $0
    exit 0
}

function child_call ()
{
    if [ "$QARRAY_CMD" == "" ] ; then
	error "in child_call $SGE_TASK_ID, no QARRAY_CMD environment given"
    fi
    if [ ! -e $QARRAY_CMD ] ; then
	error "can't locate command file $QARRAY_CMD from there: "`pwd`
    fi
    # extract line $SGE_TASK_ID $QARRAY_CMD and execute it via the SHELL defined for qsub (option -S), NOT THE USER SHELL ($SGE_O_SHELL)
    # by default SHELL is /bin/bash 
    # .. thus for qarray command file with module load requests, option -S /bin/tcsh is required if user SHELL is not bash
    # otherwise environement will be incompatible
    head -$SGE_TASK_ID $QARRAY_CMD | tail -1 | $SHELL
    exit $?
}

if [ "$SGE_TASK_ID" == "" -o "$SGE_TASK_ID" == "undefined" ] ; then 
    father_call "$@"
else
    child_call "$@"
fi
