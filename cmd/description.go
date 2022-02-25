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
	"encoding/json"
	"io/ioutil"
	"log"
	"os"
	"strings"

	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

type Organism struct {
	Vernacular string
	Binominal  string
	FunFact    string
	Taxonomy   string
}

////////////////////////////////////////////////////////////////////////////////////////////////////

type ChromlengthAssembly struct {
	Karyotype         string
	ScaffoldN50       string
	NumberOfScaffolds string
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// data struct for the decoded data
// all fields must be exportable
type Data struct {
	Organism            Organism
	ChromlengthAssembly ChromlengthAssembly
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: documentation assembly description
// descriptionCmd represents the description command
var descriptionCmd = &cobra.Command{
	Use:   "description",
	Short: "Collect assembly information.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

Parse JSON file containing description of genome assemblies.

Collect information and concatenate.
`,

	Example: ``,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	Run: func(cmd *cobra.Command, args []string) {

		// json file
		content, err := ioutil.ReadFile(inDir + "/" + species)
		if err != nil {
			log.Fatal("Error when opening file: ", err)
		}

		// unmarshall the data into features
		var features Data
		err = json.Unmarshal(content, &features)
		if err != nil {
			log.Fatal("Error during Unmarshal(): ", err)
		}

		// write extracted
		writeSpeciesFeatures(features)

	},
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	assemblyCmd.AddCommand(descriptionCmd)

	// flags

}

////////////////////////////////////////////////////////////////////////////////////////////////////

// write positions
func writeSpeciesFeatures(features Data) {

	// declare io
	f, err := os.OpenFile(outDir+"/"+outFile, os.O_APPEND|os.O_WRONLY|os.O_CREATE, 0666)

	if err != nil {
		panic(err)
	}

	defer f.Close()

	// declare writer
	w := bufio.NewWriter(f)

	// writing
	_, err = w.WriteString(
		features.Organism.Binominal + "," +
			strings.ReplaceAll(features.Organism.Vernacular, ",", "") + "," +
			features.ChromlengthAssembly.Karyotype + "," +
			strings.ReplaceAll(features.ChromlengthAssembly.ScaffoldN50, ",", "") + "," +
			strings.ReplaceAll(features.ChromlengthAssembly.NumberOfScaffolds, ",", "") +
			"\n")
	if err != nil {
		panic(err)
	}

	// flush writer
	w.Flush()
}

////////////////////////////////////////////////////////////////////////////////////////////////////
