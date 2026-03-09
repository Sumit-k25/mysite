function downloadImage(url, name) {
	var xhr = new XMLHttpRequest();
	xhr.open("GET", url, true);
	xhr.responseType = "blob";

	xhr.onload = function() {
		var urlCreator = window.URL || window.webkitURL;
		var imageUrl = urlCreator.createObjectURL(this.response);
		var tag = document.createElement('a');
		tag.href = imageUrl;
		tag.download = name+'.png';
		document.body.appendChild(tag);
		tag.click();
		document.body.removeChild(tag);
	};

	xhr.send();
}

function checkAmbiguity() {
	const sequence = document.getElementsByName("seq")[0].value.toUpperCase();
	const ambiguousCodes = ["B", "J", "O", "U", "X", "Z", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0", "!", "@", "#", "$", "%", "^", "*", "(", ")", "-", "+", "=", "[", "]", "{", "}", ";", ":", "'", "\"", "\\", "|", ",", "<", ".", ">", "/", "?", "`", "~", "&", " ", "\t", "\n"];
	const isValid = true;

	if(sequence == "") {
		alert("Input a protein sequence to proceed!");
		isValid = false;

	}else{
		for (let i = 0; i < sequence.length; i++) {
        if (ambiguousCodes.includes(sequence[i])) {
			alert("Ambiguous code detected in protein sequence: " + sequence[i]);
			isValid = false;
			break;
	        }
	    }
	}

	if (isValid){
	    document.getElementById("dataForm").submit();
		const loaderContainer = document.querySelector('.loading-screen');
		loaderContainer.style.display = 'flex';
	}
}

// Selection of all checkboxes
const check = document.getElementById("select");
const changeName = document.getElementById("select-all-functions");
check.addEventListener('change', (event) => {
	const changeName = document.getElementById("select-all-functions");

	var ele=document.getElementsByName('cal');

	if (check && check.type === "checkbox" && check.checked) {
		for(var i=0; i<ele.length; i++){
			if(ele[i].type=='checkbox')
			ele[i].checked=true;
		}
		changeName.innerText = "Deselect all options";
	}
	else{
		for(var i=0; i<ele.length; i++){
			if(ele[i].type=='checkbox')
				ele[i].checked=false;
		}
		changeName.innerText = "Select all options";
	}
});

document.addEventListener("DOMContentLoaded", function() {
    var sampleSequenceButton = document.getElementById("sampleSequenceButton");
    var seqInput = document.getElementById("seq");

    sampleSequenceButton.addEventListener("click", function() {
        var sampleSequence = "HYYLYRHVNKIT";
        seqInput.value = sampleSequence;
    });
});

// Bulk processing logic
function handleBulkUpload() {
    const fileInput = document.getElementById('seq_file');
    const form = document.getElementById('dataForm');

    if (fileInput.files.length > 0) {
        form.action = "/bulk_sequence_list";
        form.submit();
        document.querySelector('.loading-screen').style.display = 'flex';
    } else {
        alert("Please upload a .txt file!");
    }
}

// Ensure the page starts with "Deselect All" since items are pre-checked
document.addEventListener("DOMContentLoaded", function() {
    const selectAllCheck = document.getElementById("select");
    const label = document.getElementById("select-all-functions");
    
    if (selectAllCheck) {
        selectAllCheck.checked = true;
        label.innerText = "Deselect all options";
    }
});










