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
	"bufio"
	"fmt"
	"io/ioutil"
	"log"
	"os"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	sequenceID      []string
	outFastaID      string
	chromosomeField = "Linear-Chromosome"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// collectIDCmd represents the collectID command
var collectIDCmd = &cobra.Command{
	Use:   "collectID",
	Short: "Collect ID",
	Long:  `Collect ID & reformats for ` + chalk.White.Color("Chromosome List File") + `for ` + chalk.White.Color("European Nucleotide Archive") + ` submission`,

	Example: ``,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// collect ids
		collectID(fastaFile)

		// write
		writeID()
	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	fastaCmd.AddCommand(collectIDCmd)

	// flags
	collectIDCmd.Flags().StringVarP(&outFastaID, "out", "o", "", "Out file. If empty, results will be written to stdout")
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// read file & collect sequences
func collectID(readFile string) {

	// open an input file, exit on error
	contentFile, err := ioutil.ReadFile(readFile)
	if err != nil {
		log.Fatal("Error opending input file :", err)
	}

	// check whether file exists to avoid appending
	if fileExist(outFastaID) {
		os.Remove(outFastaID)
	}

	// mount data string
	dataFasta := strings.NewReader(string(contentFile))

	// fasta.Reader requires a known type template to fill
	// with FASTA data. Here we use *linear.Seq.
	template := linear.NewSeq("", nil, alphabet.DNAredundant)
	readerFasta := fasta.NewReader(dataFasta, template)

	// make a seqio.Scanner to simplify iterating over a
	// stream of data.
	scanFasta := seqio.NewScanner(readerFasta)

	// iterate through each sequence in a multifasta and
	// examine the ID, description and sequence data.
	for scanFasta.Next() {
		// get the current sequence and type assert to *linear.Seq.
		// while this is unnecessary here, it can be useful to have
		// the concrete type.
		sequence := scanFasta.Seq().(*linear.Seq)

		sequenceID = append(sequenceID, sequence.ID)
	}

	if err := scanFasta.Error(); err != nil {
		log.Fatal(err)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write positions
func writeID() {

	if outFastaID == "" {

		// printing
		for i := 0; i < len(sequenceID); i++ {
			fmt.Println(
				sequenceID[i] + "\t" +
					strings.Replace(sequenceID[i], "HiC_scaffold_", "", -1) + "\t" +
					chromosomeField,
			)
		}
	} else {

		// declare io
		f, err := os.OpenFile(outDir+"/"+outFastaID, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

		if err != nil {
			panic(err)
		}

		defer f.Close()

		// declare writer
		w := bufio.NewWriter(f)

		// writing
		for i := 0; i < len(sequenceID); i++ {
			_, err = w.WriteString(
				sequenceID[i] + "\t" +
					strings.Replace(sequenceID[i], "HiC_scaffold_", "", -1) + "\t" +
					chromosomeField + "\n",
			)

			if err != nil {
				panic(err)
			}
		}

		// flush writer
		w.Flush()
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
