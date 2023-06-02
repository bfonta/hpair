package main

import (
	"fmt"
	"os"
	"flag"
	"errors"
	"strings"
	"strconv"
	fp "path/filepath"
	ex "os/exec"
	"regexp"
)

func check(e error) {
    if e != nil {
        panic(e)
    }
}

func create_dir(d string) {
	if _, err := os.Stat(d); errors.Is(err, os.ErrNotExist) {
		os.Mkdir(d, 0700) //os.Mkdir(d, os.ModeDir)
	}
}

func run_hpair(in string, out string) {
    exe := "./run"

	fmt.Println("Command: " + exe + " " + in + " " + out)
    cmd := ex.Command(exe, in, out)

    _, err := cmd.Output()
    if err != nil {
        fmt.Println(err.Error())
        return
    }
    // fmt.Println(string(stdout))
}

func create_input_file(tri string, template []byte, out string) {
	var tri2 string;
	if strings.Contains(tri, "m") {
		tri2 = strings.Replace(tri, "m", "-", -1)
	} else {
		tri2 = tri
	}
	
	template_lines := strings.Split(string(template), "\n")
	for i, line := range template_lines {
		if strings.Contains(line, "C_3") {
			template_lines[i] = "C_3      =  " + tri2 + ".0D0";
		}
	}

	output_lines := strings.Join(template_lines, "\n")
	err := os.WriteFile(out, []byte(output_lines), 0644)
	check(err)
}

func extract_cross_sections(out string) ([]string, []string) {
	f, err := os.ReadFile(out)
	check(err)

	var matches_born [][]string;
	var matches_nlo [][]string;
	f_lines := strings.Split(string(f), "\n")
	re := regexp.MustCompile("^.+?\\(   ([0-9]+?)\\.(.+?)E(.+?)  .*?([0-9])\\.(.+)E(.+) \\) PB$")
	for _, line := range f_lines {
		if strings.Contains(line, "SIGMA_BORN") {
			//fmt.Println(line)
			matches_born = re.FindAllStringSubmatch(line, -1)
		}
		if strings.Contains(line, "SIGMA_NLO") {
			//fmt.Println(line)
			matches_nlo = re.FindAllStringSubmatch(line, -1)
		}
	}
	return matches_born[0], matches_nlo[0]
	//return []string{"TRASH", "1", "3321", "9", "2", "4322", "6"}, []string{"TRASH", "5", "21321", "9", "8", "7522", "6"}
}

func convert_to_float(arr []string) (string, string) {
	if len(arr) != 7 {
        panic("The array must have length 7! It has instead length " + strconv.Itoa(len(arr)))
    }

	var val_str string = arr[1] + "." + arr[2] + "E" + arr[3]
	var err_str string = arr[4] + "." + arr[5] + "E" + arr[6]

	return val_str, err_str
}

func save_cross_sections(out string, born_vals []string, born_errs []string, nlo_vals []string, nlo_errs []string) {
	f, err := os.Create(out)
	check(err)
    defer f.Close()

	fmt.Fprintln(f, strings.Join(born_vals, `,`) + " # Born xsection values")
	fmt.Fprintln(f, strings.Join(born_errs, `,`) + " # Born xsection errors")
	fmt.Fprintln(f, strings.Join(nlo_errs,	`,`) + " # NLO xsection values")
	fmt.Fprintln(f, strings.Join(nlo_errs,	`,`) + " # NLO xsection errors")
}

func main() {
	flag.Usage = func() {
	    fmt.Fprintf(os.Stderr, "Usage: go run create_files.go 1 2 3\n")
	}   
	flag.Parse()
	parameters := flag.Args()

	base_dir := fp.Join(string(os.Getenv("HOME")), "hpair")
	hpair_in_dir := fp.Join(base_dir, "inputs/")
	hpair_out_dir := fp.Join(base_dir, "outputs/")
	create_dir(hpair_in_dir)
	create_dir(hpair_out_dir)


	// open template file
	template, err := os.ReadFile("hpair.template.in")
	check(err)
    // defer template.Close()

	var sigma_born, error_born []string
	var sigma_nlo, error_nlo []string
	for _, tri := range parameters {

		hpair_in  := fp.Join(hpair_in_dir,  "hpair_" + tri + ".in")
		hpair_out := fp.Join(hpair_out_dir, "hpair_" + tri + ".out")
		
		create_input_file(tri, template, hpair_in)

		run_hpair(hpair_in, hpair_out)

		born, nlo := extract_cross_sections(hpair_out)

		val, error := convert_to_float(born)
		sigma_born = append(sigma_born, val)
		error_born = append(error_born, error)

		val, error = convert_to_float(nlo)
		sigma_nlo = append(sigma_nlo, val)
		error_nlo = append(error_nlo, error)
	}

	save_cross_sections(fp.Join(hpair_out_dir, "results.txt"),
		sigma_born, error_born, sigma_nlo, error_nlo)
}
