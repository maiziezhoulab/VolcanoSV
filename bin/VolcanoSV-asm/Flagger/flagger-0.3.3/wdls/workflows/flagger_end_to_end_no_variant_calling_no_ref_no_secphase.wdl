version 1.0

import "flagger.wdl" as flagger_t
import "flagger_preprocess_no_variant_calling.wdl" as preprocess_t


workflow FlaggerEndToEndNoVariantCallingNoRefNoSecphase{
    input{
        File assemblyFastaGz
        File readAlignmentBam
        Float maxReadDivergence
        String sampleName
        String suffix
        File fai
        Int threadCount = 30
    }
    call preprocess_t.runFlaggerPreprocess as preprocess{
        input:
            bam = readAlignmentBam,
            assemblyFastaGz = assemblyFastaGz,
            maxDivergence = maxReadDivergence,
            threadCount = threadCount
    }
    call flagger_t.runFlagger as flagger{
        input:
            coverageGz = preprocess.correctedCovGz,
            highMapqCoverageGz = preprocess.correctedHighMapqCovGz,
            fai = fai,
            sampleName = sampleName,
            suffix = suffix,
            covFloat = preprocess.meanCorrectedCoverageFloat,
            threadCount = threadCount
    }
    
    output {
        # flagger preprocess files
        Float meanCoverageFloat = preprocess.meanCorrectedCoverageFloat
        File covGz = preprocess.correctedCovGz
        File highMapqCovGz = preprocess.correctedHighMapqCovGz
        File excludedReadIdsText = preprocess.excludedReadIdsText

        # flagger outputs for all alignments
        File finalBed = flagger.finalBed
        File miscFilesTarGz = flagger.miscFilesTarGz
        File pdf = flagger.pdf
    }
}
