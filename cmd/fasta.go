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
	fastaFile string
)

////////////////////////////////////////////////////////////////////////////////////////////////////
// fastaCmd represents the fasta command
var fastaCmd = &cobra.Command{
	Use:   "fasta",
	Short: "Handle fasta operations.",
	Long: `Handle fasta operations, such as:

collect: collect sequence IDs.
convert: convert ` + chalk.Yellow.Color("fasta") + ` files.
segregate: segregate scaffolds into individual files.
sequence: extract sequences from assemblies.
`,

	Example: ``,
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func fastaScanner(inputFile []byte) seqio.Scanner {

	// mount data string
	dataFasta := strings.NewReader(string(inputFile))

	// fasta.Reader requires a known type template to fill with FASTA data. Here we use *linear.Seq
	template := linear.NewSeq("", nil, alphabet.DNAredundant)
	readerFasta := fasta.NewReader(dataFasta, template)

	// make a seqio.Scanner to simplify iterating over a stream of data
	scanFasta := seqio.NewScanner(readerFasta)

	// return
	return *scanFasta
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(fastaCmd)

	// persistent flags
	fastaCmd.PersistentFlags().StringVarP(&fastaFile, "fasta", "f", "", "Fasta file")

}

////////////////////////////////////////////////////////////////////////////////////////////////////
