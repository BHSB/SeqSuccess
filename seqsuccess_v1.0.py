import numpy as np
import pandas as pd
import datetime

#Open and read csv's, filter all non-robotic batch numbers from df1 and save as df3
##FUTURE: Try open using relative directory positions not full dir path.
datapath = r"C:\Users\BHSB\Dropbox\Programming\Python\Projects\MolGen\SeqSucess\data"
df1 = pd.read_csv(datapath + '\\seq17_full.csv') #'seqdata_test_5.csv'
df2 = pd.read_csv(datapath + '\\seqbatches17_full.csv') #'seqbatches_test_5.csv'

seqworkno = []
for index, row in df2.iterrows():
    if row['SeqNumbers'] not in seqworkno:
        seqworkno.append(row['SeqNumbers'])

df3 = df1[df1['Workbatch'].isin(seqworkno)] #Should be just robotic sequencing batches.

#batches with no QQ's
batches_no_qq = []
#batches that have NTC QQ but has failed
batches_failed_ntc = []
batches_passed_ntc = []
#number of batches at start of analysis
#number of batches excluded from analysis due to lack or NTC's and failed NTC's

def find_ntc(batch):
    temp_df = df1[df1['Workbatch'] == batch]

    qq_dict = {}

    for index, row in temp_df.iterrows():
        if row['SampleID'].startswith('Q'):
            qq_dict.update({index : row['SampleID']})

    if len(qq_dict) >= 1:
        return max(qq_dict.keys())
    elif len(qq_dict) == 0:
        batches_no_qq.append(batch)
        return "No NTC"

def ntc_passed(batch, index):
    if df1.iloc[index]['Result'].startswith("Fail") or df1.iloc[index]['Result'].startswith('#'):
        batches_passed_ntc.append(batch)
        return True
    else:
        batches_failed_ntc.append(batch)
        return False

def ntc_check(batch):
    qq_result = find_ntc(batch)

    if qq_result == "No NTC":
        return "No NTC"
    elif qq_result > 1:
        return ntc_passed(batch, qq_result)

for num in seqworkno:
    ntc_check(num)

print(batches_no_qq, batches_failed_ntc)

print(len(seqworkno), len(batches_no_qq), len(batches_failed_ntc), len(batches_passed_ntc))

print(len(seqworkno) - len(batches_no_qq) - len(batches_failed_ntc))

test_df = df1[df1['Workbatch'] == 1706190]

test_df.shape

find_ntc(1706190)

def num_sample(batch):
    temp_df = df1[df1['Workbatch'] == batch]
    return temp_df.shape[0]

def num_fails(batch):
    temp_df = df1[df1['Workbatch'] == batch]
    count = 0
    for index, row in temp_df.iterrows():
        if row['Result'].startswith('Fail'):
            count += 1
    return count

def num_retest(batch):
    temp_df = df1[df1['Workbatch'] == batch]
    count = 0
    for index, row in temp_df.iterrows():
        if row['Result'].startswith('#'):
            count += 1
    return count

num_sample(1706190)

def num_patients(batch):
    #Excluding positive controls (QQ's)
    temp_df = df1[df1['Workbatch'] == batch]

    patient_ids = []

    for index, row in temp_df.iterrows():
        if row['SampleID'].startswith('EX') and row['SampleID'] not in patient_ids:
            patient_ids.append(row['SampleID'])

    return patient_ids, len(patient_ids)

def num_controls(batch):
    temp_df = df1[df1['Workbatch'] == batch]

    qq_ids = []

    #if find_ntc(batch) == 1

    for index, row in temp_df.iterrows():
        if row['SampleID'].startswith('Q') and row['SampleID'] not in qq_ids:
            qq_ids.append(row['SampleID'])
    return qq_ids, len(qq_ids)

num_patients(1706190)


def unique_pts(batch):
    temp_df = df1[df1['Workbatch'] == batch]

    unique_pt = []

    for index, row in temp_df.iterrows():
        if row['SampleID'] not in unique_pt:
            unique_pt.append(row['SampleID'])
    return unique_pt

def patient_fail(batch):
    temp_df = df1[df1['Workbatch'] == batch]
    total_pts = num_patients(batch)[0] + num_controls(batch)[0]
    total_pts.pop()

    fail_count = 0
    pass_count = 0
    total_fail = 0
    patient_fail_count = 0

    for sample in total_pts:
        pt_df = temp_df[temp_df['SampleID'] == sample]
        for index, row in pt_df.iterrows():
            if row['Result'].startswith('Fail') or row['Result'].startswith('#'):
                fail_count += 1
            else:
                pass_count += 1
        if fail_count == 0:
            pass_count = 0
            fail_count = 0
        elif fail_count / (fail_count + pass_count) >= 0.5:
            patient_fail_count += 1
            total_fail += fail_count
            pass_count = 0
            fail_count = 0
        else:
            pass_count = 0
            fail_count = 0
    return total_fail, patient_fail_count

if 1700003 in batches_passed_ntc:
    print("Yes")
else:
    print("No")

def seq_success_rate(batch):

    total_samples = num_sample(batch) - 1
    total_fails = num_fails(batch) + num_retest(batch)
    deduct_fails = patient_fail(batch)[0]

    success_rate = ((total_samples - (total_fails - deduct_fails)) / total_samples) * 100

    return(total_samples, total_fails, deduct_fails, success_rate)

seq_success_rate(1706190)

#Function to create dataframes to create output csv files
def output1(): #Batches with an NTC that has passed
    print("Processing...")
    final_data = []

    for i, batch in enumerate(batches_passed_ntc):
        temp_dict = {"1. SF No." : batch, \
        "2. No. Samples" : seq_success_rate(batch)[0], \
        "3. No. Fails" : seq_success_rate(batch)[1], \
        "4. No. Patient Fails" : patient_fail(batch)[1], \
        "5. Removed Fails" : seq_success_rate(batch)[2], \
        "6. Success Rate" : seq_success_rate(batch)[3]}
        final_data.append(temp_dict)
        if i == 210:
            print("Data compiling complete")
    return final_data

def output2(): #Batches with a failed NTC or no NTC
    print("Processing...")
    final_data = []

    for batch in batches_no_qq:
        temp_dict = {"1. SF No." : batch, "2. No. Samples" : num_sample(batch),\
        "3. No Fails" : num_fails(batch) + num_retest(batch)}
        final_data.append(temp_dict)

    for batch in batches_failed_ntc:
        temp_dict = {"1. SF No." : batch, "2. No. Samples" : num_sample(batch),\
        "3. No Fails" : num_fails(batch) + num_retest(batch)}
        final_data.append(temp_dict)

    print("Data compiling complete.")
    return final_data

def total_avg():
    avg_total = 0
    for index, row in output1.iterrows():
        avg_total += row['6. Success Rate']
    print(output1.shape[0])
    return (avg_total / output1.shape[0])

output1 = pd.DataFrame(output1())

output2 = pd.DataFrame(output2())

date_now = datetime.datetime.now()
datestamp = str(date_now.year) + str(date_now.month) + str(date_now.day)

output1.to_csv("C:\\Users\\BHSB\\Desktop\\Seq files\\" + "output1_" + datestamp + ".csv")

output2.to_csv("C:\\Users\\BHSB\\Desktop\\Seq files\\" + "output2_" + datestamp + ".csv")
