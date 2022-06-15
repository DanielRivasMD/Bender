/*
Copyright © 2021 Daniel Rivas <danielrivasmd@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
*/
package cmd

import (
	"bufio"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	identity float64
	length   float64
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// lociCmd represents the loci command
var lociCmd = &cobra.Command{
	Use:   "loci",
	Short: "Identify loci from similarity search result.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

Identify & filter loci form similarity search result, based on identity & length.
`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(κ *cobra.Command, args []string) {

		// execute logic
		genomicPositionsCollect(inDir + "/" + species)

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(lociCmd)

	// flags
	lociCmd.Flags().Float64VarP(&identity, "identity", "i", 80, "Percentage of sequence identity")
	lociCmd.Flags().Float64VarP(&length, "length", "l", 400, "Length of sequence alignment")

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// read file, filter records & write
func genomicPositionsCollect(readFile string) {

	// open an input file, exit on error
	inputFile, ε := os.Open(readFile)
	if ε != nil {
		log.Fatal("Error opening input file : ", ε)
	}

	// check whether file exists to avoid appending
	if fileExist(outDir + "/" + species) {
		os.Remove(outDir + "/" + species)
	}

	// headers := []string{
	// 	"qseqid",
	// 	"sseqid",
	// 	"pident",
	// 	"length",
	// 	"mismatch",
	// 	"gapopen",
	// 	"qstart",
	// 	"qend",
	// 	"sstart",
	// 	"send",
	// 	"evalue",
	// 	"bitscore",
	// }

	// scanner.Scan() advances to the next token returning false if an error was encountered
	scanner := bufio.NewScanner(inputFile)

	for scanner.Scan() {

		// tab separated records
		records := strings.Split(scanner.Text(), "\t")

		// collect patterns
		pIdent, _ := strconv.ParseFloat(records[2], 64)
		alignLen, _ := strconv.ParseFloat(records[3], 64)

		// filter criteria
		if pIdent >= identity && alignLen >= length {
			// write
			writeGenomicPositions(outDir+"/"+species, records)
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write positions
func writeGenomicPositions(outFile string, records []string) {

	// declare io
	ƒ, ε := os.OpenFile(outFile, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

	if ε != nil {
		panic(ε)
	}

	defer ƒ.Close()

	// declare writer
	ϖ := bufio.NewWriter(ƒ)

	// writing
	_, ε = ϖ.WriteString(
		records[0] + "\t" +
			records[1] + "\t" +
			records[2] + "\t" +
			records[3] + "\t" +
			records[4] + "\t" +
			records[5] + "\t" +
			records[6] + "\t" +
			records[7] + "\t" +
			records[8] + "\t" +
			records[9] + "\t" +
			records[10] + "\t" +
			records[11] +
			"\n")

	if ε != nil {
		panic(ε)
	}

	// flush writer
	ϖ.Flush()
}

////////////////////////////////////////////////////////////////////////////////////////////////////
