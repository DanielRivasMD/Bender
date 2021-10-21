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
	"io/ioutil"
	"log"
	"math"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/biogo/biogo/alphabet"
	"github.com/biogo/biogo/io/seqio"
	"github.com/biogo/biogo/io/seqio/fasta"
	"github.com/biogo/biogo/seq/linear"

	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	syncytin   identified // identified struct
	scaffoldID string
	startCoor  float64
	endCoor    float64
	hood       float64
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// syncytin features
type identified struct {
	scaffoldIdent string
	positionIdent position
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// positions
type position struct {
	startPos float64
	endPos   float64
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// SequenceCmd represents the Sequence command
var SequenceCmd = &cobra.Command{
	Use:   "Sequence",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	Example: ``,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// scaffold
		syncytin.scaffoldIdent = scaffoldID

		// positions
		syncytin.positionIdent.minMax(startCoor, endCoor)

		// declare file output
		if outFile == "" {
			outFile = sequenceOut(species)
		}

		// execute logic
		collectCoordinates(inDir + "/" + species)

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	AssemblyCmd.AddCommand(SequenceCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	SequenceCmd.Flags().StringVarP(&scaffoldID, "scaffold", "S", "", "Scaffold ID")
	SequenceCmd.Flags().Float64VarP(&startCoor, "start", "b", 0., "Start coordinate")
	SequenceCmd.Flags().Float64VarP(&endCoor, "end", "e", 0., "End coordinate")
	SequenceCmd.Flags().Float64VarP(&hood, "hood", "n", 20000, "Neighborhood to look into")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// pass struct as reference to update
func (position *position) updateMinMax(num1, num2 float64) {
	position.startPos = math.Min(num1, num2)
	position.endPos = math.Max(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// determine minimum / maximum
func (position *position) minMax(num1, num2 float64) {
	position.updateMinMax(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// define output file
func sequenceOut(readFile string) string {

	outFile = readFile
	reg := regexp.MustCompile(`HiC*`)
	res := reg.FindStringIndex(outFile)
	outFile = outFile[0:res[0]]

	// assembly directory
	outFile = outFile +
		syncytin.scaffoldIdent + "_" +
		strconv.FormatFloat(syncytin.positionIdent.startPos, 'f', 0, 64) + "_" +
		strconv.FormatFloat(syncytin.positionIdent.endPos, 'f', 0, 64) +
		".fasta"
	return outFile
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// read file & collect sequences
func collectCoordinates(readFile string) {

	// open an input file, exit on error
	contentFile, err := ioutil.ReadFile(readFile)
	if err != nil {
		log.Fatal("Error opending input file :", err)
	}

	// check whether file exists to avoid appending
	if fileExist(outFile) {
		os.Remove(outFile)
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

		// find scaffold
		if sequence.ID == syncytin.scaffoldIdent {

			// cast coordinates
			startI := int64(syncytin.positionIdent.startPos - hood)
			endI := int64(syncytin.positionIdent.endPos + hood)

			id := syncytin.scaffoldIdent + "_" + strconv.FormatFloat(syncytin.positionIdent.startPos, 'f', 0, 64) + "_" + strconv.FormatFloat(syncytin.positionIdent.endPos, 'f', 0, 64)
			// find coordinates
			targatSeq := linear.NewSeq(id, sequence.Seq[startI:endI], alphabet.DNA)

			// write candidate
			writeFasta(outFile, targatSeq)
		}

	}

	if err := scanFasta.Error(); err != nil {
		log.Fatal(err)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write positions
func writeFasta(outFile string, sequence *linear.Seq) {

	// declare io
	f, err := os.OpenFile(outFile, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

	if err != nil {
		panic(err)
	}

	defer f.Close()

	// declare writer
	w := fasta.NewWriter(f, 10000)

	// writing
	_, err = w.Write(sequence)

	if err != nil {
		panic(err)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
