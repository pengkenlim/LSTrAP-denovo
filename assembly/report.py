#setting sys.path for importing modules
import os
import sys
if __name__ == "__main__":
        abspath= os.getcwd()
        parent_module= os.path.join(abspath.split("HSS-Trans")[0], "HSS-Trans")
        sys.path.insert(0, parent_module)
        
import json
import pandas as pd
from numpy import round

def generate_from_json_log(logpath, reportpath,step=1):
    '''As the function name suggests, 
    this function will create a HTML report by 
    extracting information from json logfile in the output directory'''
    #opening logfile and loading in dictionary
    with open(logpath , "r") as f:
        log_content= json.load(f)
    
    #extracting info from log.json's dictionary and assigning them as variables
    
    #information relavent to Step_1
    consensus_v_cds_data= [['Consensus theshold', 'Number of CDS']] + [[idx+1, n_cds] for idx, n_cds in enumerate([stats[0] for stats in log_content["Step_1"]["consensus"]["stats"].values()])]
    taxid = log_content["Step_1"]["run_info"].get("taxid")
    sci_name = log_content["Step_1"]["run_info"].get("sci_name")
    n_total_acc = log_content["Step_1"]["run_info"].get("n_total_acc")
    prelim_command= log_content["Step_1"]["run_info"].get("command_issued")
    prelim_start_time= log_content["Step_1"]["run_info"].get("init_time")
    processed_table= {"Accession processed" :log_content["Step_1"].get("processed_acc").keys(), 
    "Number of CDS assembled": log_content["Step_1"].get("processed_acc").values()}
    processed_table= pd.DataFrame.from_dict(processed_table)
    
    CT_table= {"Consensus threshold": list(log_content["Step_1"]["consensus"]["stats"].keys()),
    "Number of CDS" : [stats[0] for stats in log_content["Step_1"]["consensus"]["stats"].values()], 
    "Average CDS length": [stats[1] for stats in log_content["Step_1"]["consensus"]["stats"].values()],
    "GC content (%)": [stats[2] for stats in log_content["Step_1"]["consensus"]["stats"].values()]}
    CT_table= pd.DataFrame.from_dict(CT_table)
    optimal_CT= log_content["Step_1"]["consensus"].get("optimal")
    
    optimal_CT_string=f"<p>User is reccomended to use <b>CT{optimal_CT}</b> as preliminary assembly for step 2.</p>"
    
    if step == 1:
        Consensus_threshold_for_preliminary_assembly="NA"
        cluster_command="NA"
        cluster_start_time="NA"
        ps_cutoff ,n_failed , n_total  ="NA" , "NA" , "NA"
        histo_data="NA"
        sc_data="NA"
        optimal_k , sc_max ="NA" , "NA"
        med_cluster_size, min_cluster_size , max_cluster_size = "NA", "NA", "NA"
        cluster_size_data= "NA"
    elif step == 2:
        #information relavent to cluster_accession.py
        Consensus_threshold_for_preliminary_assembly =  log_content["Step_2"]["run_info"].get("Consensus_threshold_for_preliminary_assembly")
        cluster_command = log_content["Step_2"]["run_info"].get("command_issued")
        cluster_start_time = log_content["Step_2"]["run_info"].get("init_time")
        ps_cutoff = round(log_content["Step_2"]["qc"].get("threshold"), 2)
        n_failed = len(log_content["Step_2"]["qc"].get("failed"))
        n_total = len(log_content["Step_2"]["qc"].get("total"))
    
        processed_acc = log_content["Step_2"].get("processed_acc")
        histo_data = [["Accession",'Pseudoalignment rate (%)'] ]+ [[key,val] for key, val in processed_acc.items() if type(val)==float ]
    
        k_sc = log_content["Step_2"]["kmeans"].get("s_coeficient")
        sc_data = [["k","Silhouette coefficient"]] + [[int(key), val] for key, val in k_sc.items()]
    
        optimal_k= log_content["Step_2"]["kmeans"].get("cluster_assignment_stats")[0]
        sc_max = round(log_content["Step_2"]["kmeans"].get("cluster_assignment_stats")[1], 2)
        med_cluster_size=log_content["Step_2"]["kmeans"].get("cluster_assignment_stats")[2]
        mean_cluster_size = log_content["Step_2"]["kmeans"].get("cluster_assignment_stats")[3]
        min_cluster_size = log_content["Step_2"]["kmeans"].get("cluster_assignment_stats")[4]
        max_cluster_size = log_content["Step_2"]["kmeans"].get("cluster_assignment_stats")[5]

        cluster_size_data = [["Cluster", "Number of accessions", { "role": "style" }]]
        for i in range(0,optimal_k):
            if i%2 == 0:
                t_string= "#ff9912"
            else:
                t_string= "#4682b4"
            cluster_size_data+= [[f"Cluster {str(i)}", len(log_content["cluster"]["kmeans"]["cluster_assignment_dict"].get(str(i))),t_string ]]
        cluster_size_data += [["QC failed", n_failed , "#bebebe"]]

    
    html_string= f'''
            <html>
	<head>
		<title>HSS-Trans report</title>
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
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
	<script type="text/javascript">
      google.charts.load("current", {{packages:["corechart"]}});
      google.charts.setOnLoadCallback(drawChart);
      function drawChart() {{
        var data = google.visualization.arrayToDataTable({histo_data}
		  );

        var options = {{
          title: 'Pseudoalignment of Accessions Against Preliminary Assembly',
          legend: {{ position: 'none' }},
		  histogram: {{ bucketSize: 5}},
		  vAxis: {{title: "Number of accessions"}},
		  hAxis: {{title: "Pseudoalignment rate (%)"}}
        }};

        var chart = new google.visualization.Histogram(document.getElementById('PS_HISTO'));
        chart.draw(data, options);
      }}
    </script>
    <script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
    <script type="text/javascript">
      google.charts.load('current', {{'packages':['corechart']}});
      google.charts.setOnLoadCallback(drawChart);

      function drawChart() {{
        var data = google.visualization.arrayToDataTable({sc_data});

       var options = {{
          title: 'Clustering Performance of K-means Clustering Iterations at Different ks',
          legend: {{ position: 'none' }},
		  vAxis: {{title: "Silhouette coefficient"}},
		  hAxis: {{title: "Number of clusters (k)"}}
        }};

        var chart = new google.visualization.LineChart(document.getElementById('SC_CURVE'));

        chart.draw(data, options);
      }}
    </script>
    	<script type="text/javascript" src="https://www.gstatic.com/charts/loader.js"></script>
	<script type="text/javascript">
		google.charts.load("current", {{packages:["corechart"]}});
		google.charts.setOnLoadCallback(drawChart);
		function drawChart() {{
		  var data = google.visualization.arrayToDataTable({cluster_size_data});

		  var options = {{
			title: "Size of clusters",
			legend: {{position: "none" }},
		  }};
		  var chart = new google.visualization.BarChart(document.getElementById("CLUSTER_SIZE_BAR"));
		  chart.draw(data, options);
	  }}
	</script>
	</head>
	<body>
		<h1 style="background-color:#D6EEEE;">HSS-Trans Report</h1>
		<h2><u>H</u>igh-throughput <u>S</u>ample <u>S</u>election pipeline for <u>Trans</u>criptome assembly</h2>
		<p>For more information, please visit our <a href="https://github.com/pengkenlim/HSS-Trans">GitHub repository</a>.</p>
		<h2 style="background-color:#D6EEEE;" >Step 1. Assembling Draft CDSs (reduced but high-confidence assembly)</h2>
		<h3>Run Info</h3>
		<table>
			<tr>
				<th style = "text-align: right;" >Organism NCBI TaxID:</th>
				<td>{taxid}</td>
			</tr>
			<tr>
				<th style = "text-align: right;" >Scientific name:</th>
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
        {CT_table.to_html(index=False)}
		<h2 style="background-color:#D6EEEE;">Step 2. Selection of representative accessions for transcriptome assembly</h2>
        <h3>Run Info</h3>
		<table>
			<tr>
				<th style = "text-align: right;" >Command issued :</th>
				<td>{cluster_command}</td>
			</tr>
            <tr>
				<th style = "text-align: right;" >Preliminary assembly consensus threshold:</th>
				<td>{Consensus_threshold_for_preliminary_assembly}</td>
			</tr>
			<tr>
				<th style = "text-align: right;" >Initial start time:</th>
				<td>{cluster_start_time}</td>
			</tr>
		</table>
        <h3>Accession Quality Control</h3>
        <div id="PS_HISTO" style="width: 900px; height: 500px;"></div>
        <p>Pseudoalignment threshold of <b>{ps_cutoff}%</b> was used for quality control (QC). <b>{n_failed}</b> out of {n_total} accessions failed QC and was excluded out of pipeline.</p>
        <h3>Clustering of Accessions using K-means Algorithm</h3>
        <div id="SC_CURVE" style="width: 900px; height: 500px;"></div>
        <p>Optimal K-means iteration determined to be at <b>k={optimal_k}</b> with a silhouette coefficient of <b>{sc_max}</b>.</p>
        <div id="CLUSTER_SIZE_BAR" style="width: 900px; height: 900px;"></div>
		<h2 style="background-color:#D6EEEE;">Step 3. Modular concatenation of multi-kmer and transcriptome-profile-specific assemblies.</h2>	
		<h3>Format to be determined</h3>
	</body>
</html>'''
    #creating report file
    with open(reportpath, "w") as f:
        f.write(html_string)
