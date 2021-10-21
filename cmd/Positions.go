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
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	identity float64
	length   float64
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: documentation assembly positions
// PositionsCmd represents the Positions command
var PositionsCmd = &cobra.Command{
	Use:   "Positions",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// execute logic
		genomicPositionsCollect(inDir + "/" + species)

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	AssemblyCmd.AddCommand(PositionsCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	PositionsCmd.Flags().Float64VarP(&identity, "identity", "i", 80, "Percentage of sequence identity")
	PositionsCmd.Flags().Float64VarP(&length, "length", "l", 400, "Length of sequence alignment")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// read file, filter records & write
func genomicPositionsCollect(readFile string) {

	// open an input file, exit on error
	inputFile, readErr := os.Open(readFile)
	if readErr != nil {
		log.Fatal("Error opening input file : ", readErr)
	}

	// check whether file exists to avoid appending
	if fileExist(outFile) {
		os.Remove(outFile)
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
	f, err := os.OpenFile(outFile, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

	if err != nil {
		panic(err)
	}

	defer f.Close()

	// declare writer
	w := bufio.NewWriter(f)

	// writing
	_, err = w.WriteString(
		records[0] + " " +
			records[1] + " " +
			records[2] + " " +
			records[6] + " " +
			records[7] + " " +
			records[10] +
			"\n")
	if err != nil {
		panic(err)
	}

	// flush writer
	w.Flush()
}

////////////////////////////////////////////////////////////////////////////////////////////////////
