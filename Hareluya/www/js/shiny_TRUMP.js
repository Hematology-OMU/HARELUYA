$(function() {
  
  const log = console.log;

  // const

  const ID_TARGET = 'distPlot';
  const CHART_TITLE = 'Calculated OS by Random Survival Forest';
  const TABLE_TITLE = 'Probability of OS at 1-year after allo-HSCT';
  const LEGEND_TITLE = 'Allo-HSCT procedures';
  const LOADING_IMG_PATH = 'img/loading.gif';
  const TOOLTIP_OFFSET = 20;

  const LOADING_TEXT = 
	'Calculating the probability of OS for {total} allo-HSCT procedures.<br>' +
	'The number of remaining calculation is {rest} allo-HSCT procedures.'

  // member

  $('#' + ID_TARGET)
    .append($('<div>').attr('id', 'chart_panel')
      .append($('<div>').attr('id', 'chart')
        .append($('<p>').attr('class', 'title').text(CHART_TITLE))
        .append($('<div>').attr('id', 'chart_image'))
      )
      .append($('<div>').attr('id', 'loading')
        .append($('<div>').attr('id', 'loading_box')
          .append($('<p>').attr('id', 'loading_text'))
          .append($('<img>').attr('src', 'img/loading.gif'))
        )
      )
    );

  $('#' + 'chart-ranking')
    .append($('<p>').attr('class', 'title').text(TABLE_TITLE))
    .append($('<dl>').attr('id', 'values_table')
  );

  $('#' + "chart-legend")
    .append($('<p>').attr('class', 'title').text(LEGEND_TITLE))
    .append($('<ul>').attr('id', 'chart_series')
      .append($('<li>').attr('id', 'R-MRD'))
      .append($('<li>').attr('id', 'M-MRD'))
      .append($('<li>').attr('id', 'R-MUD'))
      .append($('<li>').attr('id', 'M-MUD'))
      .append($('<li>').attr('id', 'R-UCB'))
      .append($('<li>').attr('id', 'M-UCB'))
      .append($('<li>').attr('id', 'R-Haplo-CY'))
      .append($('<li>').attr('id', 'M-Haplo-CY'))
      .append($('<li>').attr('id', 'R-Haplo-nonCY'))
      .append($('<li>').attr('id', 'M-Haplo-nonCY'))
    );

  /**
   * inetractive-graph
   */
  const chart = c3.generate({
    bindto: '#chart_image',
    data: {
      x: 'x',
      columns: [],
      colors: {
        'R-MRD': '#e44d93',
        'R-MUD': '#00a0de',
        'R-UCB': '#ee7b1a',
        'R-Haplo-CY': '#a9cc51',
        'R-Haplo-nonCY': '#9b7cb6',
        'M-MRD': '#e5171f',
        'M-MUD': '#0078ba',
        'M-UCB': '#814721',
        'M-Haplo-CY': '#019a66',
        'M-Haplo-nonCY': '#522886'
      },
      done: onLoad
    },
    point: {
      show: false
    },
    axis: {
      font: '12px sans-serif',
      x: {
        min: 0,
        max: 25,
        label: {
          text: 'Months after transplantation',
          position: 'outer-center'
        },
        tick: {
          values: [0, 5, 10, 15, 20]
        }
      },
      y: {
        min: 0.0,
        max: 1.0,
        label: {
          text: 'Overall Survival',
          position: 'outer-middle'
        },
        tick: {
          values: [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
        }
      }
    },
    padding: {
      bottom: 20,
    },
    legend: {
      show: false
    },
    tooltip: {
      format: {
        title: function (x, index) {
          return floatFixed(x, 2) + ' Month';
        },
        name: function (name, ratio, id, index) {
          return name;
        },
        value: function (value, ratio, id) {
          return floatFixed(value, 3);
        }
      },
      position: function (data, width, height, element) {
        const mouse = d3.mouse(element);
        const chartX = document.querySelector('#chart_image').getBoundingClientRect().x;
        const chartWidth = document.querySelector('#chart_image').getBoundingClientRect().width;
        const xgridFocusX = document.querySelector('.c3-xgrid-focus').getBoundingClientRect().x - chartX;
        const tooltipWidth = document.querySelector('.c3-tooltip-container').getBoundingClientRect().width;
        const x =
          ((xgridFocusX + TOOLTIP_OFFSET + tooltipWidth) <= chartWidth) ?
          (xgridFocusX + TOOLTIP_OFFSET) :
          (xgridFocusX - (TOOLTIP_OFFSET + tooltipWidth));
        const y = mouse[1];
        return {top: y, left: x}
      }
    }
  });

  /**
   * display-loading
  */

  // helper

  function floatFixed(value, num) {
    const digit = 10**num;
    return ((Math.round(value * digit)) / digit).toFixed(num);
  }

  function loadingText(m) {
    return LOADING_TEXT
      .replace('{total}', String(m.total))
      .replace('{rest}', String(m.rest));
  }
  
  // action
  
  function onLoad() {
    {
      const xs = Object.entries(chart.x())[0][1];
      const index = xs.findIndex(function (e) {
        return (e > 12);
      }) - 1;
      const es = Object.entries(chart.data()).map(function ([_, v]) {
        return {
          name: v.id,
          value: v.values[index].value
        };
      });
      const ess = es.sort(function (lhs, rhs) {
        return rhs.value - lhs.value;
      });

      $('#values_table').empty();
      ess.forEach(function (e, i) {
        $('#values_table')
          .append($('<dt>').text(i + 1))
          .append($('<dd>').text(e.name))
          .append($('<dd>').text(floatFixed(e.value, 3)))
      });
    }

    {
      $('#chart_series li').empty();
      Object.entries(chart.data.colors()).map(function ([name, color]) {
        $('#' + name)
          .append($('<span>')
            .css('backgroundColor', color)
            .on('click', toggleSeries)
          )
          .append($('<span>')
            .text(name)
            .on('click', toggleSeries)
          );
      });
    }
  }
  
  function showLoading(m) {
    console.log(m);
    $('#chart').fadeOut('normal', function () {
      $('#loading').fadeIn('normal');
    });
    $('#loading_text').html(loadingText(m));
  }

  function toggleSeries() {
    const e = $(this).parent();
    const id = e.attr('id');
    if (e.hasClass('not-show')) {
      chart.show(id);
    } else {
      chart.hide(id);
    }
    e.toggleClass('not-show');
  }

  function updateLoading(m) {
    console.log(m);
    $('#loading_text').html(loadingText(m));
  }

  function hideLoading(m) {
    console.dir(m);
    m = JSON.stringify(m);
    m = JSON.parse(m, function (k, v) {
      if (typeof v === 'string') {
        if (v.match(/^[-]?(\d*)(\.\d+)?$/)) {
          v = parseFloat(v);
        }
      }
      return v;
    });
    chart.load({
        columns: m.data,
        done: onLoad
    });
    $('#loading').fadeOut('normal', function () {
      $('#chart').fadeIn('normal');
    });
  }

  // handler

  function calcHandler(m) {
    if (m.status === 'start') {
      showLoading(m);
    } else if (m.status === 'processing') {
      updateLoading(m); 
    } else if (m.status === 'end') {
      hideLoading(m); 
    }
  }

  // bind

  Shiny.addCustomMessageHandler('calc', calcHandler);

  // init

  chart.load({
    url: 'data/data.csv',
    done: onLoad
  });
  
});
