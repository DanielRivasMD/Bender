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
	"os"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"
	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// segregateCmd represents the segregate command
var segregateCmd = &cobra.Command{
	Use:   "segregate",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// segregate ids
		segregateID(inDir + "/" + fastaFile)

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
	contentFile, ε := ioutil.ReadFile(readFile)
	if ε != nil {
		log.Fatal("Error opending input file :", ε)
	}

	// check whether file exists to avoid appending
	if fileExist(outFastaID) {
		os.Remove(outFastaID)
	}

	// mount data string
	dataFasta := strings.NewReader(string(contentFile))

	// fasta.Reader requires a known type template to fill with FASTA data. Here we use *linear.Seq
	template := linear.NewSeq("", nil, alphabet.DNAredundant)
	readerFasta := fasta.NewReader(dataFasta, template)

	// make a seqio.Scanner to simplify iterating over a stream of data
	scanFasta := seqio.NewScanner(readerFasta)

	// iterate through each sequence in a multifasta and examine the ID, description and sequence data
	for scanFasta.Next() {
		// get the current sequence and type assert to *linear.Seq while this is unnecessary here, it can be useful to have the concrete type
		sequence := scanFasta.Seq().(*linear.Seq)

		fmt.Printf(">%s\n%s\n", sequence.ID, sequence.Seq)

	}

	if ε := scanFasta.Error(); ε != nil {
		log.Fatal(ε)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
