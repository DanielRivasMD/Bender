/*
Copyright © 2020 Daniel Rivas <danielrivasmd@gmail.com>

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
	"os"
	"strings"

	"github.com/spf13/cobra"
	"github.com/spf13/pflag"
	"github.com/spf13/viper"
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var (
	cfgPath     string
	cfgFile     string
	verboseBool string
)

////////////////////////////////////////////////////////////////////////////////////////////////////

var rootCmd = &cobra.Command{
	Use:   "bender",
	Short: "A robot to handle Slurm Genomic jobs",
	Long: `Daniel Rivas <danielrivasmd@gmail.com>

"Good news everyone!"
Bender is a robot for automation on
Genomic jobs in Slurm systems.
"It's highly addictive!"

Bender creates a convinient command line interphase
with built-in and accessible documentation`,
	Version: "v0.2",
	Example: `
Bender help`,

	////////////////////////////////////////////////////////////////////////////////////////////////////

	PersistentPreRun: func(cmd *cobra.Command, arg []string) {

		initializeConfig(cmd, cfgPath, strings.TrimSuffix(cfgFile, ".toml"))

	},

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func Execute() {
	if err := rootCmd.Execute(); err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////

func init() {

	////////////////////////////////////////////////////////////////////////////////////////////////////

	// TODO: standarize flag order on call
	// persistent flags
	rootCmd.PersistentFlags().StringVarP(&cfgPath, "configPath", "p", ".", "Path to config file")

	rootCmd.PersistentFlags().StringVarP(&cfgFile, "configFile", "c", "", "Config file")

	// rootCmd.PersistentFlags().StringP("inDir", "i", "", "Directory where input files are located")
	// viper.BindPFlag("inDir", rootCmd.PersistentFlags().Lookup("inDir"))

	// rootCmd.PersistentFlags().StringP("outDir", "o", "", "Output directory. Creates if not exitst")
	// viper.BindPFlag("outDir", rootCmd.PersistentFlags().Lookup("outDir"))

	rootCmd.PersistentFlags().StringVarP(&verboseBool, "verbose", "v", "false", "Verbosity switch")

	////////////////////////////////////////////////////////////////////////////////////////////////////

}

////////////////////////////////////////////////////////////////////////////////////////////////////

func initializeConfig(cmd *cobra.Command, configPath string, configName string) error {

	// initialize viper
	v := viper.New()

	// collect config path & file from persistent flags
	v.AddConfigPath(configPath)
	v.SetConfigName(configName)

	// read the config file
	if err := v.ReadInConfig(); err != nil {
		// okay if there isn't a config file
		if _, ok := err.(viper.ConfigFileNotFoundError); !ok {
			// return an error if we cannot parse the config file
			return err
		}
	}

	// bind flags to viper
	bindFlags(cmd, v)

	return nil
}

////////////////////////////////////////////////////////////////////////////////////////////////////

// bind each cobra flag to its associated viper configuration
func bindFlags(cmd *cobra.Command, v *viper.Viper) {

	cmd.Flags().VisitAll(func(f *pflag.Flag) {

		// apply the viper config value to the flag when the flag is not set and viper has a value
		if !f.Changed && v.IsSet(f.Name) {
			val := v.Get(f.Name)
			cmd.Flags().Set(f.Name, fmt.Sprintf("%v", val))
		}
	})
}

////////////////////////////////////////////////////////////////////////////////////////////////////
