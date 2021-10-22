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
	"github.com/spf13/cobra"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	assembly string
	species  string
	outFile  string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// TODO: documentation assembly
// assemblyCmd represents the assembly command
var assemblyCmd = &cobra.Command{
	Use:   "assembly",
	Short: "Handle assembly operations.",
	Long: `Handle assembly operations, such as:

description: parse assembly features.
loci: read diamond assembly search output & filter results.
search: perform similarity search.
sequence: extract sequences from assemblies.
synteny: .
`,

	Example: ``,

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {

	////////////////////////////////////////////////////////////////////////////////////////////////////

	rootCmd.AddCommand(assemblyCmd)

	// persistent flags
	assemblyCmd.PersistentFlags().StringVarP(&assembly, "assembly", "a", "", "Assembly file")
	assemblyCmd.PersistentFlags().StringVarP(&species, "species", "s", "", "Species file")
	assemblyCmd.PersistentFlags().StringVarP(&outFile, "outfile", "o", "", "Out file. If empty it will be defined by input")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////
