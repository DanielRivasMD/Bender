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
	fileOut  string     // infered from input
	syncytin identified // identified struct
	scaffold string
	start    float64
	end      float64
	suffix   string
	hood     int64
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// syncytin features
type identified struct {
	scaffold  string
	positions position
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// positions
type position struct {
	start float64
	end   float64
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
		syncytin.scaffold = scaffold

		// positions
		syncytin.positions.parseMinMax(start, end)

		// declare file output
		fileOut = defineOut(assembly, suffix)

		// execute logic
		collectCoordinates(inDir + "/" + "DNAzoo" + "/" + assembly)

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// define output file
func defineOut(readFile, suffix string) string {

	var path string
	if suffix == "" {
		path = "candidate"
	} else {
		path = "insertion"
	}

	fileOut = readFile
	reg := regexp.MustCompile(`HiC*`)
	res := reg.FindStringIndex(fileOut)
	fileOut = fileOut[0:res[0]]

	fileOut = inDir + "/" + path + "/" +
		fileOut +
		syncytin.scaffold + "_" +
		strconv.FormatFloat(syncytin.positions.start, 'f', 0, 64) + "_" +
		strconv.FormatFloat(syncytin.positions.end, 'f', 0, 64) +
		suffix + ".fasta"
	return fileOut
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// pass struct as reference to update
func (position *position) minMax(num1, num2 float64) {
	position.start = math.Min(num1, num2)
	position.end = math.Max(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// parse values & determine minimum / maximum
func (position *position) parseMinMax(num1, num2 float64) {
	position.minMax(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// read file & collect sequences
func collectCoordinates(readFile string) {

	contentFile, err := ioutil.ReadFile(readFile) // the file is inside the local directory
	if err != nil {
		log.Fatal("Error opending input file :", err)
	}

	// check whether file exists to avoid appending
	if fileExist(fileOut) {
		os.Remove(fileOut)
	}

	// mount data string
	dataFasta := strings.NewReader(string(contentFile))

	// fasta.Reader requires a known type template to fill
	// with FASTA data. Here we use *linear.Seq.
	template := linear.NewSeq("", nil, alphabet.DNAredundant)
	readerFasta := fasta.NewReader(dataFasta, template)

	// Make a seqio.Scanner to simplify iterating over a
	// stream of data.
	scanFasta := seqio.NewScanner(readerFasta)

	// Iterate through each sequence in a multifasta and examine the
	// ID, description and sequence data.
	for scanFasta.Next() {
		// Get the current sequence and type assert to *linear.Seq.
		// While this is unnecessary here, it can be useful to have
		// the concrete type.
		sequence := scanFasta.Seq().(*linear.Seq)

		// find scaffold
		if sequence.ID == syncytin.scaffold {

			// cast coordinates
			startI := int64(syncytin.positions.start)
			endI := int64(syncytin.positions.end)

			// extract neighborhood
			switch suffix {
			case "_upstream":
				endI = startI
				startI = startI - hood
			case "_downstream":
				startI = endI
				endI = endI + hood
			}

			id := syncytin.scaffold + "_" + strconv.FormatFloat(syncytin.positions.start, 'f', 0, 64) + "_" + strconv.FormatFloat(syncytin.positions.end, 'f', 0, 64)
			// find coordinates
			targatSeq := linear.NewSeq(id, sequence.Seq[startI:endI], alphabet.DNA)

			// write candidate
			writeFasta(fileOut, targatSeq)
		}

	}

	if err := scanFasta.Error(); err != nil {
		log.Fatal(err)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write positions
func writeFasta(fileOut string, sequence *linear.Seq) {

	// declare io
	f, err := os.OpenFile(fileOut, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

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

// fileExist checks if a file exists and is not a directory before
// try using it to prevent further errors
func fileExist(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	AssemblyCmd.AddCommand(SequenceCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags
	SequenceCmd.Flags().StringVarP(&scaffold, "scaffold", "i", "", "Scaffold ID")
	SequenceCmd.Flags().Float64P("start", "s", 0., "Start coordinate")
	SequenceCmd.Flags().Float64P("end", "e", 0., "End coordinate")
	SequenceCmd.Flags().StringVarP(&suffix, "suffix", "x", "", "Suffix")
	SequenceCmd.Flags().Int64P("hood", "n", 10000, "Neighborhood to look into")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////
