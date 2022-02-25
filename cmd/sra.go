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
	"github.com/ttacon/chalk"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

// declarations
var ()

////////////////////////////////////////////////////////////////////////////////////////////////////

// sraCmd represents the sra command
var sraCmd = &cobra.Command{
	Use:   "sra",
	Short: "Handles SRA operations.",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

Collect files from Short Read Archive (SRA)
and check the state of the downloads.

Reconstruct binary SRA files to fastq.
`,

	Example: `
` + chalk.Cyan.Color("bender") + ` sra help`,
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {
	rootCmd.AddCommand(sraCmd)

	// flags

}

////////////////////////////////////////////////////////////////////////////////////////////////////
