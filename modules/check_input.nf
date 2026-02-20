def check_input(genomes_tsv_path) {

    new File(genomes_tsv_path).withReader { reader ->
        def header = reader.readLine()
        def columns = header.split("\t")
        def expectedColumns = ["genome", "species", "representative", "genome_is_representative", "fasta_path"]

        // Check that columns are expected
        if (columns.size() != expectedColumns.size()) {
            throw new IllegalStateException("Error: File at ${genomes_tsv_path} \
            does not have exactly ${expectedColumns.size()} columns")
        }
        if (!columns.equals(expectedColumns)) {
            throw new IllegalStateException("Error: Column headers do not \
            match expected headers: ${expectedColumns.join(", ")}")
        }

        def speciesMap = [:].withDefault { [representativeCount: 0, totalCount: 0] }
        def genomeSet = new HashSet()

        reader.eachLine { line ->
            def fields = line.split("\t")
            if (fields.size() != expectedColumns.size()) {
                throw new IllegalStateException("Error: Inconsistent number of columns in line: ${line}")
            }
            
            def genome = fields[0]
            def species = fields[1]
            def genomeIsRepresentative = fields[3] as Integer
            def fastaPath = new File(fields[4])

            // Check that all genome names are unique
            if (!genomeSet.add(genome)) {
                throw new IllegalStateException("Error: Duplicate genome ID found: ${genome}")
            }
            
            // Check that genome name is POSIX compliant
            if (genome.contains("/") || genome.contains("\0")) {
                throw new IllegalStateException("Error: Genome name contains illegal characters: ${genome}")
            }

            // Count number of representative genomes per species
            speciesMap[species].totalCount++
            if (genomeIsRepresentative == 1) {
                speciesMap[species].representativeCount++
            }

            // Check that genome fasta file exists and is not empty
            if (!fastaPath.exists()) {
                throw new IllegalArgumentException("File not found for genome: ${genome} at ${fastaPath}")
            }        
            if (fastaPath.length() == 0) {
                throw new IllegalArgumentException("File is empty: ${fastaPath}")
            }
        }
        
    }
}

