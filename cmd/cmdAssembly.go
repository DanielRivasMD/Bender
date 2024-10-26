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
	"math"
	"strconv"

	"github.com/spf13/cobra"
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var (
	assembly   string
	species    string
	library    string
	libraryDir string
	outFile    string
	syncytin   identified // identified struct
	scaffold   string
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

// annotations
type annotation struct {
	scaffoldIdent   string
	classIdent      string
	positionIdent   position
	scoreInt        int64
	strandSense     string
	attributeStruct attribute
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// attributes
type attribute struct {
	ID     string
	Name   string
	Alias  string
	Parent string
	Target string
	Note   string
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// assemblyCmd represents the assembly command
var assemblyCmd = &cobra.Command{
	Use:   "assembly",
	Short: "Handle assembly operations.",
	Long: chalk.Green.Color("Daniel Rivas <danielrivasmd@gmail.com>") + `

Handle assembly operations, such as:

	` + chalk.Magenta.Color("database") + `:	build database.
	` + chalk.Magenta.Color("description") + `:	parse assembly features.
	` + chalk.Magenta.Color("filter") + `:		filter assembly list.
	` + chalk.Magenta.Color("loci") + `:		read diamond assembly search output & filter results.
	` + chalk.Magenta.Color("search") + `:		perform similarity search.
	` + chalk.Magenta.Color("synteny") + `:	collect syntenic information.
`,

	Example: ``,
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

// parse values & determine minimum / maximum
func (position *position) parseMinMax(str1, str2 string) {
	num1, _ := strconv.ParseFloat(str1, 64)
	num2, _ := strconv.ParseFloat(str2, 64)

	position.minMax(num1, num2)
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(assemblyCmd)

	// persistent flags
	assemblyCmd.PersistentFlags().StringVarP(&assembly, "assembly", "a", "", "Assembly file")
	assemblyCmd.PersistentFlags().StringVarP(&species, "species", "p", "", "Species file")
	assemblyCmd.PersistentFlags().StringVarP(&outFile, "outfile", "o", "", "Out file. If empty it will be defined by input")

	assemblyCmd.PersistentFlags().StringVarP(&library, "library", "l", "", "Library to search against")
	assemblyCmd.PersistentFlags().StringVarP(&libraryDir, "libraryDir", "L", "", "Library directory")

	assemblyCmd.PersistentFlags().StringVarP(&scaffold, "scaffold", "S", "", "Scaffold")
	assemblyCmd.PersistentFlags().Float64VarP(&startCoor, "start", "s", 0, "Start coordinate")
	assemblyCmd.PersistentFlags().Float64VarP(&endCoor, "end", "e", 0, "End coordinate")
	assemblyCmd.PersistentFlags().Float64VarP(&hood, "hood", "n", 0, "Neighborhood to look into")

}

////////////////////////////////////////////////////////////////////////////////////////////////////
