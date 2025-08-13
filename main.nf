#!/usr/bin/env nextflow

// import the workflow code from the hello.nf file
include { CutRun } from './workflows/CutRun.nf'



// declare input parameter
//params.metadata = '/Users/stubbe/Desktop/metadata.csv'

workflow {

//Calls the cleanRead workflow
    CutRun()
}