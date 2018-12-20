# Heuristics file in use at the NYU's Center for Brain Imaging
# (Our data come from a Siemens Prisma 3T scanner).
# Author: Pablo Velasco
# Date: 03/28/2018

import os

# TO DO: What do we do about pairs of scans "original" + "normalize" (when normalizing
#        for the coil sensitivity profile, but also keeping the original?)
#
#        Maybe I should use ".lower()" ???
#
#        How to treat Mag and Phase images for the same run?
#        Right now, it treats them as two different series, so they have
#        a different run number.  However, they really are the same run no.
#        I can check the "image_type" and see if it has 'P' in it (phase), and
#        if so, go to the previous one with 'M' and use the same run number?

def create_key(template, outtype=('nii.gz','dicom'), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')
    return (template, outtype, annotation_classes)


def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where

    allowed template fields - follow python string module:

    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    # anat:
    t1_scout = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-Scout_run-{item:02d}_T1w')
    t1 = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T1w')
    t2 = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_T2w')
    pd_bias_body = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-biasBody_run-{item:02d}_PD')
    pd_bias_receive = create_key('{bids_subject_session_dir}/anat/{bids_subject_session_prefix}_acq-biasReceive_run-{item:02d}_PD')
    # func:
    # For the functional runs, we want to have a different key for each task.  Since we don't know a-priori what the task
    #  names will be, we don't create the different keys here, but we'll create them dynamically, as needed.
    # diffusion:
    dwi = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_dwi')
    dwi_sbref = create_key('{bids_subject_session_dir}/dwi/{bids_subject_session_prefix}_acq-{acq}_run-{item:02d}_sbref')
    # fmaps:
    fmap_topup = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_dir-{direction}_run-{item:02d}_epi')    # for "topup"
    fmap_topup_AP = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_dir-AP_run-{item:02d}_epi')    # for "topup", AP
    fmap_topup_PA = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_dir-PA_run-{item:02d}_epi')    # for "topup", PA
    fmap_topup_RL = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_dir-RL_run-{item:02d}_epi')    # for "topup", RL
    fmap_topup_LR = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-fMRI_dir-LR_run-{item:02d}_epi')    # for "topup", LR
    fmap_gre_mag = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-GRE_run-{item:02d}_magnitude')       # GRE fmap
    fmap_gre_phase = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-GRE_run-{item:02d}_phasediff')     #
    fmap_dwi = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-{direction}_run-{item:02d}_epi')
    fmap_dwi_AP = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-AP_run-{item:02d}_epi')
    fmap_dwi_AP_sbref = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-AP_run-{item:02d}_sbref')
    fmap_dwi_PA = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-PA_run-{item:02d}_epi')
    fmap_dwi_RL = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-RL_run-{item:02d}_epi')
    fmap_dwi_LR = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_acq-dwi_dir-LR_run-{item:02d}_epi')
    #phoenix_doc = create_key('{bids_subject_session_dir}/fmap/{bids_subject_session_prefix}_phoenix')

    info = {t1_scout:[], t1:[], t2:[], pd_bias_body:[], pd_bias_receive:[],
            dwi:[], dwi_sbref:[],
            fmap_topup:[], fmap_topup_AP:[], fmap_topup_PA:[], fmap_topup_RL:[], fmap_topup_LR:[],
            fmap_gre_mag:[], fmap_gre_phase:[],
            fmap_dwi:[], fmap_dwi_AP:[], fmap_dwi_PA:[], fmap_dwi_RL:[], fmap_dwi_LR:[], fmap_dwi_AP_sbref:[]}
            #phoenix_doc:[]}

    for idx, s in enumerate(seqinfo):
        #pdb.set_trace()
        # s is a namedtuple with fields equal to the names of the columns
        # found in the dicominfo.tsv file

        acq = ''

        ###   T1w   ###
        # 1) Auto-Align Head scout (the original images, not the derived):
        #    look for "AA*Scout" or "AA*scout" in protocol_name:
        if ('scout' in s.protocol_name.lower()) and (s.is_derived == False):
            info[t1_scout].append({'item': s.series_id})
        # 2) High resolution T1w:
        # single volume, protocol name including T1, T1w, MPRAGE, MP-RAGE, MPR,...
        if ((s.dim4 == 1) and ( ('t1' in s.protocol_name.lower()) or
                                (('mp' in s.protocol_name.lower()) and ('rage' in s.protocol_name.lower())) or
                                ('mpr' in s.protocol_name.lower()) ) and
                              ('fl' in s.sequence_name)):
            # check the PE ('_PA' or '_rev' means 'reversed'):
            if ('_AP' in s.protocol_name):
                acq = 'highresAP'
            elif ('_PA' in s.protocol_name):
                acq = 'highresPA'
            elif ('_rev' in s.protocol_name):
                acq = 'highresrev'
            else:
                acq = 'highres'
            info[t1].append({'item': s.series_id, 'acq': acq})
        ###   T2w   ###
        # single volume, protocol name including T2, T2w, TSE, SPACE, SPC:
        if (s.dim4 == 1) and (('T2' in s.protocol_name) or
                              ('TSE' in s.protocol_name) or
                              ('SPACE' in s.protocol_name) or
                              ('SPC' in s.protocol_name) ):
            # check the PE ('_PA' or '_rev' means 'reversed'):
            if ('_AP' in s.protocol_name):
                acq = 'highresAP'
            elif ('_PA' in s.protocol_name):
                acq = 'highresPA'
            elif ('_rev' in s.protocol_name):
                acq = 'highresrev'
            else:
                acq = 'highres'
            info[t2].append({'item': s.series_id, 'acq': acq})

        ###   PD   ###
        # BIAS images  (for coil sensitivity estimation) are typically PD-weighted
        if ('bias' in s.protocol_name.lower()) and ('tfl3d' in s.sequence_name):
            if ('body' in s.protocol_name.lower()):
                info[pd_bias_body].append({'item': s.series_id})
            else:
                info[pd_bias_receive].append({'item': s.series_id})

        ###   FUNCTIONAL   ###
        # We want to make sure the _SBRef and phase series (if present) are labeled
        #  the same as the main (magnitude) image.  So we only focus on the magnitude
        #  series (to exclude phase images) and more than 3 volumes (to exclude _SBRef)
        #  and then we search if the phase and/or _SBRef are present.
        if ((s.dim4 >= 4) and ('epfid2d' in s.sequence_name)
                          and ('M' in s.image_type)
                          and (s.series_description[-6:].lower() != '_sbref')):

            ###   functional -- check PE direction   ###
            # ('_PA' or '_rev' means 'reversed')
            # (This is common to all functional runs, so we check it first)
            if ('_AP' in s.protocol_name):
                acq = 'AP'
            elif ('_PA' in s.protocol_name):
                acq = 'PA'
            elif ('_rev' in s.protocol_name):
                acq = 'rev'
            else:
                acq = 'normal'


            ###   functional -- check task name   ###
            task = ''
            if ('rest' in s.protocol_name.lower()):
                task = 'rest'
            elif ('face' in s.protocol_name.lower()):
                task = 'face'
            elif ('conflict' in s.protocol_name.lower()):
                task = 'conflict'
            elif ('gambl' in s.protocol_name.lower()):
                task = 'gamble'
            elif ('task' in s.protocol_name.lower()):
                # we want to capture what comes after "task", until the next
                #    dash ("-") or underscore ("_"):
                task = s.protocol_name.lower().split('task')[1]
                # remove any initial "-" or "_":
                if ((task[0] == '-') or (task[0] == '_')):
                    task = task[1:]
                # discard anything after the next "-" or "_":
                task = task.split('_')[0]
                task = task.split('-')[0]
                # remove any spaces we might have:
                task = task.replace(" ", "")
            else:
                task = 'TASK'    # fallback.  BIDS requires an alphanumeric string (no spaces...)

            # Below, for both magnitude and phase cases, we check if a template is defined
            #    with that specific name. That way we can assure the run number (run-01,
            #    run-02, ...) only apply for that specific task, not for just any
            #    task. E.g., if I have two tasks, with two different runs for each
            #    task name they should be called: task-foo_run-01, task-foo_run-02,
            #    task-boo_run-01,...    as opposed to task-boo_run-03 because there
            #    were 2 previous "task" runs.

            ###   functional -- is phase image present?   ###
            # At least for Siemens systems, if magnitude/phase was selected, the
            #  phase images come as a separate series immediatelly following the
            #  magnitude series.
            # (note: make sure you don't check beyond the number of elements...)

            if (idx+1 < len(seqinfo)) and ('P' in seqinfo[idx+1].image_type):
                # we have a magnitude/phase pair:

                # dictionary keys specific for this task type:
                mykey_mag = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_rec-magnitude_run-{item:02d}_bold" % task)
                mykey_pha = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_rec-phase_run-{item:02d}_bold" % task)
                try:
                    # check if "info" already has this key by trying to append to it.
                    info[mykey_mag].append({'item': s.series_id,                'acq': acq})
                    info[mykey_pha].append({'item': seqinfo[idx + 1].series_id, 'acq': acq})
                except KeyError:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey_mag] = [{'item': s.series_id,                'acq': acq}]
                    info[mykey_pha] = [{'item': seqinfo[idx + 1].series_id, 'acq': acq}]

            else:
                # we only have a magnitude image

                # dictionary key specific for this task type:
                mykey = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_run-{item:02d}_bold" % task)
                try:
                    # check if "info" already has this key by trying to append to it.
                    info[mykey].append({'item': s.series_id, 'acq': acq})
                except KeyError:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey] = [{ 'item': s.series_id, 'acq': acq }]

            ###   SB REF   ###
            # here, within the functional run code, check to see if the
            #  previous run protocol name ended in _SBREF, to assign the
            #  same task name and --if possible-- same run number.
            if (idx > 0) and ('_sbref' in seqinfo[idx - 1].series_description.lower()):
                mykey_sb = create_key("{bids_subject_session_dir}/func/{bids_subject_session_prefix}_task-%s_acq-{acq}_run-{item:02d}_sbref" % task)
                try:
                    # check if "info" already has this key by trying to append to it.
                    info[mykey_sb].append({'item': seqinfo[idx - 1].series_id, 'acq': acq})
                except KeyError:
                    # if it doesn't, add this key, specifying the first item:
                    info[mykey_sb] = [{'item': seqinfo[idx - 1].series_id, 'acq': acq}]


        ###   FIELD MAPS   ###
        # TO-DO: how do we make sure each FM has the same number as the corresponding functional?

        # A) Spin-Echo distortion maps
        #    to be used with topup:
        #pdb.set_trace()
        if ((s.dim4 <= 3) and ('epse2d' in s.sequence_name)
                          and (  ('dist'  in s.protocol_name.lower())
                              or ('map'   in s.protocol_name.lower())
                              or ('field' in s.protocol_name.lower()) )):

            if (s.series_description[-4:] != '_SBRef'):    # sbref from MB diffusion have epse2d in
                                                           #  sequence_name, so don't include them
                ###   SE distortion: -- check PE direction   ###
                if ('_AP' in s.protocol_name):
                    info[fmap_topup_AP].append({'item': s.series_id})
                elif ('_PA' in s.protocol_name):
                    info[fmap_topup_PA].append({'item': s.series_id})
                elif ('_RL' in s.protocol_name):
                    info[fmap_topup_RL].append({'item': s.series_id})
                elif ('_LR' in s.protocol_name):
                    info[fmap_topup_LR].append({'item': s.series_id})
                else:
                    if ('_rev' in s.protocol_name):
                        direction = 'rev'
                    else:
                        direction = 'normal'

                    info[fmap_topup].append({'item': s.series_id, 'direction': direction})


            # TO-DO: do more classification here (task, etc.)



        # B) GRE fmap:
        if ('fm2d' in s.sequence_name) and (('fieldmap' in s.protocol_name.lower()) or
                                                ('field_map' in s.protocol_name.lower())):
#            pdb.set_trace()
            if (s.image_type[2] == 'M'):
                # magnitude image
                # For now, we don't need to separate according to TEs, because
                #   dcm2niix seems to create 2 NIfTI files appending "1" and "2"
                #   at the end, which is what BIDS specifies.
                # If this changes, we would have to do something, like s.TE
                info[fmap_gre_mag].append({'item': s.series_id})
            elif (s.image_type[2] == 'P'):
                # phase image:
                info[fmap_gre_phase].append({'item': s.series_id})
        # "topup" fmap:
        #
        #    dirtype = s.protocol_name.split('_')[-1]
        #    info[fmap_topup].append({'item': s.series_id, 'dir': dirtype})
#            acq = s.protocol_name.split('FieldMap_')[1] + 'AP'
#            info[dwi].append({'item': s.series_id, 'acq': acq})

        ###   DIFFUSION   ###

        # We could also check: (s.image_type[2] == 'DIFFUSION')
        if ('ep_b' in s.sequence_name):       # Siemens product diffusion sequence:

            # This is not very rigorous, but in general, diffusion runs will
            #   have more than a couple of volumes.  I'll use 5, to be safe:
            if ( s.dim4 >= 5 ):
                # this is a standard diffusion acquisition
                acq = str(s.dim4)+'vols'
                info[dwi].append({'item': s.series_id, 'acq': acq})

                # check to see if the previous run is a SBREF:
                if ( (idx > 0) and
                         (seqinfo[idx - 1].series_description[-4:] == '_SBRef') and
                         ('epse2d' in seqinfo[idx - 1].sequence_name) ):
                    info[dwi_sbref].append({'item': s.series_id, 'acq': acq})

            else:
                # this is a fmap for diffusion.

                # TO-DO: for now, just ignore the _sbref image (if present, it would
                #        be the previous run (seqinfo[idx - 1].series_description[-4:] == '_SBRef')
                #        because BIDS doesn't allow them.

                # see if we can get the orientation, for topup:
                if ('_AP' in s.protocol_name):
                    info[fmap_dwi_AP].append({'item': s.series_id})
                    if ( (idx > 0) and
                            (seqinfo[idx - 1].series_description[-4:] == '_SBRef') ):
                        info[fmap_dwi_AP_sbref].append({'item': seqinfo[idx - 1].series_id})
                elif ('_PA' in s.protocol_name):
                    info[fmap_dwi_PA].append({'item': s.series_id})
                elif ('_RL' in s.protocol_name):
                    info[fmap_dwi_RL].append({'item': s.series_id})
                elif ('_LR' in s.protocol_name):
                    info[fmap_dwi_LR].append({'item': s.series_id})
                else:
                    if ('_rev' in s.protocol_name):
                        direction = 'rev'
                    else:
                        direction = 'normal'

                    info[fmap_dwi].append({'item': s.series_id, 'direction': direction})

        ###   PHOENIX FILE   ###

        #if ('PhoenixZIPReport' in s.series_description) and (s.image_type[3] == 'CSA REPORT'):       # 
        #    # Phoenix Report:
        #    info[phoenix_doc].append({'item': s.series_id})

        #if (s.dim4 >= 99) and (('dMRI_dir98_AP' in s.protocol_name) or ('dMRI_dir99_AP' in s.protocol_name)):
        #    acq = s.protocol_name.split('dMRI_')[1].split('_')[0] + 'AP'
        #    info[dwi].append({'item': s.series_id, 'acq': acq})
        #if (s.dim4 >= 99) and (('dMRI_dir98_PA' in s.protocol_name) or ('dMRI_dir99_PA' in s.protocol_name)):
        #    acq = s.protocol_name.split('dMRI_')[1].split('_')[0] + 'PA'
        #    info[dwi].append({'item': s.series_id, 'acq': acq})
        #if (s.dim4 == 1) and (('dMRI_dir98_AP' in s.protocol_name) or ('dMRI_dir99_AP' in s.protocol_name)):
        #    acq = s.protocol_name.split('dMRI_')[1].split('_')[0]
        #    info[fmap_dwi].append({'item': s.series_id, 'dir': 'AP', 'acq': acq})
        #if (s.dim4 == 1) and (('dMRI_dir98_PA' in s.protocol_name) or ('dMRI_dir99_PA' in s.protocol_name)):
        #    acq = s.protocol_name.split('dMRI_')[1].split('_')[0]
        #   info[fmap_dwi].append({'item': s.series_id, 'dir': 'PA', 'acq': acq})

    #pdb.set_trace()
    return info