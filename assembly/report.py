#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("LSTrAP-denovo")[0], "LSTrAP-denovo")
        sys.path.insert(0, parent_module)
        
import json
import pandas as pd

def generate_from_json_log(logpath, reportpath):
    '''As the function name suggests, 
    this function will create a HTML report by 
    extracting information from json logfile in the output directory'''
    #opening logfile and loading in dictionary
    with open(logpath , "r") as f:
        log_content= json.load(f)
    
    #extracting info from log.json's dictionary and assigning them as variables
    
    #information relavent to prelim
    consensus_v_cds_data= [['Consensus theshold', 'Number of CDS']] + [[CT, n_cds] for CT, n_cds in log_content["prelim"]["consensus"]["CDS"].items()]
    taxid = log_content["prelim"]["run_info"].get("taxid")
    sci_name = log_content["prelim"]["run_info"].get("sci_name")
    n_total_acc = log_content["prelim"]["run_info"].get("n_total_acc")
    prelim_command= log_content["prelim"]["run_info"].get("command_issued")
    prelim_start_time= log_content["prelim"]["run_info"].get("init_time")
    processed_table= {"Accession processed" :log_content["prelim"].get("processed_acc").keys(), 
    "Number of CDS assembled": log_content["prelim"].get("processed_acc").values()}
    processed_table= pd.DataFrame.from_dict(processed_table)
    optimal_CT= log_content["prelim"]["consensus"]["stats"].get("CT")
    
    if log_content["prelim"]["run_var"]["consensus_threshold"]==0:
        optimal_CT_string=f"<p>An optimal CT of <b>{optimal_CT}</b> has been determined automatically.</p>"
    else:
        optimal_CT_string=f"<p>CT of <b>{optimal_CT}</b> has been defined by user.</p>"
   
    n_cds= log_content["prelim"]["consensus"]["stats"].get("n_CDS")
    avg_cds_len= log_content["prelim"]["consensus"]["stats"].get("CDS_len")
    GC = log_content["prelim"]["consensus"]["stats"].get("GC")
    
    html_string= f'''
            <html>
	<head>
		<title>LSTrAP-denovo report</title>
		<style>
			table, td, th {{
			  border: 1px solid black;
			}}

			table {{
			  border-collapse: collapse;
			  width: 50%;
			}}

			td {{
			  text-align: center;
			}}
			th {{
			  text-align: center;
			}}
			tr:nth-child(even) {{
				background-color: #dcdcdc;
			}}
		</style>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
	<script type="text/javascript">
      google.charts.load('current', {{'packages':['corechart']}});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {{
        var data = google.visualization.arrayToDataTable({consensus_v_cds_data});

        var options = {{
          title: 'Number of CDS extracted with increasing consensus thresholds',
          legend: {{ position: 'none' }},
		  vAxis: {{title: "Number of CDS"}},
		  hAxis: {{title: "Consensus Threshold (CT)", ticks: {list(range(1, len(consensus_v_cds_data)))}}}
        }};

        var chart = new google.visualization.LineChart(document.getElementById('curve_chart'));

        chart.draw(data, options);
      }}
    </script>

	</head>
	<body>
		<h1 style="background-color:#D6EEEE;">LSTrAP-<i>denovo</i> Report</h1>
		<h2><u>L</u>arge <u>S</u>cale (<i>de novo</i>) <u>T</u>ranscriptome <u>A</u>ssembly <u>P</u>ipeline using publicly available RNA-seq data.</h2>
		<p>For more information, please visit our <a href="https://github.com/pengkenlim/LSTrAP-denovo">GitHub repository</a>.</p>
		<h2 style="background-color:#D6EEEE;" >Step 1. Generating a Preliminary Assembly.</h2>
		<h3>Run Info</h3>
		<table>
			<tr>
				<th style = "text-align: right;" >Organism NCBI TaxID:</th>
				<td>{taxid}</td>
			</tr>
			<tr>
				<th style = style = "text-align: right;" >Scientific name:</th>
				<td>{sci_name}</td>
			</tr>
			<tr>
				<th style = "text-align: right;" >Total number of RNA-seq accessions:</th>
				<td>{n_total_acc}</td>
			</tr>
			<tr>
				<th style = "text-align: right;" >Command issued :</th>
				<td>{prelim_command}</td>
			</tr>
			<tr>
				<th style = "text-align: right;" >Initial start time:</th>
				<td>{prelim_start_time}</td>
			</tr>
		</table>
		<h3>Download and Independent Assembly of Single-accessions</h3>
        {processed_table.to_html(index=False)}
		<h3>Combining Single-accessions Assemblies and Extracting Consensus Coding Sequences (CDS)</h3>
		<div id="curve_chart" style="width: 900px; height: 500px"></div>
		{optimal_CT_string}
			<table style= "width: 25%;">
			<tr>
				<th colspan="2", style = "text-align: center;">Statistics of Preliminary Assembly</th>
			</tr>
			<tr>
				<th style ="text-align: right;" >Consensus threshold used to extract CDS:</th>
				<td>{optimal_CT}</th>	
			</tr>
			<tr>
				<th>Number of CDS:</th>
				<td style= " text-align: right; padding-left: 30px; padding-right: 30px;">{n_cds}</td>
			</tr>
			<tr>
				<th style ="text-align: right;" >Average CDS length (nt):</th>
				<td>{avg_cds_len}</td>
			</tr>
			<tr>
				<th style ="text-align: right;" >GC content (%):</th>
				<td>{GC}</td>
			</tr>
			</table>
		<h2 style="background-color:#D6EEEE;">Step 2. Large-scale Download, Quality Control and Transcriptome-profile-based Clustering of Accessions.</h2>
		<h3>Format to be determined</h3>
		<h2 style="background-color:#D6EEEE;">Step 3. Modular concatenation of multi-kmer and transcriptome-profile-specific assemblies.</h2>	
		<h3>Format to be determined</h3>
	</body>
</html>'''
    #creating report file
    with open(reportpath, "w") as f:
        f.write(html_string)
