/*
Copyright © 2022 Daniel Rivas <danielrivasmd@gmail.com>

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
	"fmt"
	"io/ioutil"
	"log"
	"strings"

	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// segregateCmd represents the segregate command
var segregateCmd = &cobra.Command{
	Use:   "segregate",
	Short: "Segregate scaffolds into individual files.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

Segregate scaffolds into individual files.
`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// segregate ids
		segregateID(strings.Replace(fastaFile, ".gz", "", -1))

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	fastaCmd.AddCommand(segregateCmd)

	// flags

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// read file & collect sequences
func segregateID(readFile string) {

	// open an input file, exit on error
	contentFile, ε := ioutil.ReadFile(inDir + "/" + readFile)
	if ε != nil {
		log.Fatal("Error opending input file :", ε)
	}

	// scan fasta
	scanFasta := fastaScanner(contentFile)

	// iterate through each sequence in a multifasta and examine the ID, description and sequence data
	fmt.Println("Segregating scaffolds...")
	for scanFasta.Next() {

		// get the current sequence and type assert to *linear.Seq while this is unnecessary here, it can be useful to have the concrete type
		sequence := scanFasta.Seq().(*linear.Seq)

		// assign outfile dynamically
		switch {
		case strings.Contains(readFile, "fna"):
			outFile = outDir + "/" + strings.Replace(readFile, ".fna", "", -1) + "_" + sequence.ID + ".fasta"
		case strings.Contains(readFile, "faa"):
			outFile = outDir + "/" + strings.Replace(readFile, ".faa", "", -1) + "_" + sequence.ID + ".fasta"
		case strings.Contains(readFile, "fasta"):
			outFile = outDir + "/" + strings.Replace(readFile, ".fasta", "", -1) + "_" + sequence.ID + ".fasta"
		default:
			log.Fatal("fasta format not recognized")
		}

		// write sequence
		writeFasta(sequence)
	}

	if ε := scanFasta.Error(); ε != nil {
		log.Fatal(ε)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
