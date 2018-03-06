BIOFORMATS = {
    "DATAMATRIX_EXPRESSION": {
        value: 'DATAMATRIX_EXPRESSION',
        text: 'Data matrix expression',
        validator: 'ExpressionValidator',
        hint: 'Data matrix expression is a Tab-separated values file. This file has two columns if there is only one sample, and more than two columns if there are many samples. First line is a header and must contain the names of the samples. The first column correspond to gene/probe/protein Ids from (Ensembl gene, HGNC symbol, Entrez id, Affy HG U133A probeset, Affy HG U133B probeset, Affy HG U133-PLUS-2 probeset and Affy HTA 2.0), the next columns correspond to gene expression values in numeric format from each sample.',
        url: 'http://hipathia.babelomics.org/doc/doku.php?id=expression_matrix_file_format'
    },
    "VARIANT": {
        value: 'VARIANT',
        text: 'Variant (VCF)',
        validator: 'VCFValidator',
        hint: 'VCF is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome. The format also has the ability to contain genotype information on samples for each position.',
        url: 'https://samtools.github.io/hts-specs/VCFv4.2.pdf'
    },
    "EXPERIMENTAL_DESIGN": {
        value: 'EXPERIMENTAL_DESIGN',
        text: 'Experimental design',
        validator: 'ExperimentalDesignValidator',
        hint: 'Experimental design is Tab-separated values file. This file has two columns, the first is the sample name and the second is the phenotype.',
        url: 'http://hipathia.babelomics.org/doc/doku.php?id=experimental_design_file_format'
    }
};
