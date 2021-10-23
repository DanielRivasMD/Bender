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
	"math"
	"os"
	"reflect"
	"strconv"
	"strings"

	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////
// synteny range
const (
	header = "seqid" + "," + "start" + "," + "end" + "," + "source" + "," + "score" + "," + "strand" + "," + "distance" + "," + "att_id" + "," + "att_alias" + "," + "att_note" + "," + "att_target" + "\n"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	headg bool = true // head gene
	headr bool = true // head repeatmasker

	// command line arguments
	// inDir       string = os.Args[1]
	// species       string = os.Args[2]
	// scaffoldID string = os.Args[3]
	// startCoor   string = os.Args[4]
	// endCoor     string = os.Args[5]
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// syntenyCmd represents the synteny command
var syntenyCmd = &cobra.Command{
	Use:   "synteny",
	Short: "A brief description of your command",
	Long: `A longer description that spans multiple lines and likely contains examples
and usage of using your command. For example:

Cobra is a CLI library for Go that empowers applications.
This application is a tool to generate the needed files
to quickly create a Cobra application.`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// scaffold
		syncytin.scaffoldIdent = scaffoldID

		// positions
		syncytin.positionIdent.minMax(startCoor, endCoor)

		// declare file output
		if outFile == "" {
			outFile = coordinateOut(species)
		}

		// execute logic
		annotate(inDir + "/" + species)

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(syntenyCmd)

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// flags

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// print annotations
func (annotations *annotation) print() string {
	return annotations.scaffoldIdent + "," +
		strconv.FormatFloat(annotations.positionIdent.startPos, 'f', 0, 64) + "," +
		strconv.FormatFloat(annotations.positionIdent.endPos, 'f', 0, 64) + "," +
		annotations.classIdent + "," +
		strconv.FormatInt(annotations.scoreInt, 10) + "," +
		annotations.strandSense + "," +
		strconv.FormatFloat(syncytin.positionIdent.distanceCandidate(annotations.positionIdent), 'f', 0, 64) + "," +
		annotations.attributeStruct.ID + "," +
		annotations.attributeStruct.Alias + "," +
		annotations.attributeStruct.Note + "," +
		annotations.attributeStruct.Target + "\n"

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// calculate distance from candidate
func (candidate *position) distanceCandidate(annotation position) float64 {
	out := 0.
	upstream := annotation.endPos - candidate.startPos
	downstream := annotation.startPos - candidate.endPos
	minimum := math.Min(math.Abs(upstream), math.Abs(downstream))
	if math.Abs(upstream) == math.Abs(downstream) {
		out = upstream
	} else if minimum == math.Abs(upstream) {
		out = upstream
	} else if minimum == math.Abs(downstream) {
		out = downstream
	}
	return out
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// read file & collect annotations
func annotate(species string) {

	// open an input file, exit on error
	inputFile, readErr := os.Open(species)
	if readErr != nil {
		log.Fatal("Error opending input file :", readErr)
	}

	// gene check whether file exists to avoid appending
	if fileExist(outFile + "_gene.csv") {
		os.Remove(outFile + "_gene.csv")
	}

	// repeatmasker check whether file exists to avoid appending
	if fileExist(outFile + "_repm.csv") {
		os.Remove(outFile + "_repm.csv")
	}

	// scanner.Scan() advances to the next token returning false if an error was encountered
	scanner := bufio.NewScanner(inputFile)

	for scanner.Scan() {

		// split line records by tab
		records := strings.Split(scanner.Text(), "\t")

		// collect patterns. internal values are redeclared every iteration
		annotationCollect(records)

	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// annotations := []string{
// 	"seqid",
// 	"source",
// 	"type",
// 	"start",
// 	"end",
// 	"score",
// 	"strand",
// 	"phase",
// 	"attributes",
// }

////////////////////////////////////////////////////////////////////////////////////////////////////

// collect annotations
func annotationCollect(records []string) {

	if len(records) > 1 {

		// declare annotation struct
		var annotations annotation

		// scaffold
		annotations.scaffoldIdent = records[0]

		// class / type
		annotations.classIdent = records[2]

		// positions
		annotations.positionIdent.parseMinMax(records[3], records[4])

		// score
		annotations.scoreInt, _ = strconv.ParseInt(records[5], 10, 64)

		// strand
		annotations.strandSense = records[6]

		// raw attributes
		rawAttributes := records[8]

		// segregate attributes
		attributeSegregate(rawAttributes, &annotations.attributeStruct)

		if annotations.scaffoldIdent == syncytin.scaffoldIdent &&
			annotations.positionIdent.startPos > (syncytin.positionIdent.startPos-hood) &&
			annotations.positionIdent.endPos < (syncytin.positionIdent.endPos+hood) {

			if annotations.classIdent == "gene" {

				// write
				writeSyntenyGenes(outFile, "gene", annotations)

			} else if records[1] == "repeatmasker" && annotations.classIdent == "match" {

				// write
				writeSyntenyGenes(outFile, "repm", annotations)

			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// attributes := []strings{
// 	"ID",
// 	"Name",
// 	"Alias",
// 	"Parent",
// 	"Target",
// 	"Gap",
// 	"Derives_from",
// 	"Note",
// 	"Dbxref",
// 	"Ontology_term",
// 	"Is_circular",
// }

////////////////////////////////////////////////////////////////////////////////////////////////////

// pass struct as reference to update
func attributeSegregate(rawAttributes string, attributes *attribute) {

	// collect attribute struct field names
	fields := reflect.TypeOf(*attributes)

	// collect attributes
	arrayAttributes := strings.Split(rawAttributes, ";")

	// loop over attribute string array
	for ix := 0; ix < len(arrayAttributes); ix++ {

		num := fields.NumField()
		// loop over attribute struct fields
		for i := 0; i < num; i++ {
			field := fields.Field(i)
			attributes.AddAttribute(arrayAttributes[ix], field.Name)
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// add collected attribute
func (attributes *attribute) AddAttribute(ats, field string) {
	if strings.Contains(ats, field) {
		out := strings.TrimPrefix(ats, field+"=")
		final := reflect.ValueOf(attributes).Elem()
		final.FieldByName(field).Set(reflect.ValueOf(out))
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write positions
func writeSyntenyGenes(outFile, suffixOut string, annotations annotation) {

	// declare io
	f, err := os.OpenFile(outFile+"_"+suffixOut+".csv", os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

	if err != nil {
		panic(err)
	}

	defer f.Close()

	// declare writer
	w := bufio.NewWriter(f)

	// gene header
	if headg && suffixOut == "gene" {
		headg = false
		_, err = w.WriteString(header)

		if err != nil {
			panic(err)
		}
	}

	// repeatmasker header
	if headr && suffixOut == "repm" {
		headr = false
		_, err = w.WriteString(header)

		if err != nil {
			panic(err)
		}
	}

	// writing
	_, err = w.WriteString(annotations.print())

	if err != nil {
		panic(err)
	}

	// flush writer
	w.Flush()
}

////////////////////////////////////////////////////////////////////////////////////////////////////
